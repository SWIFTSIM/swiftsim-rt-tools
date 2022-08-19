#!/usr/bin/env python3

###############################################################################
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


import swiftsimio
import matplotlib as mpl

mpl.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
import sys
import stromgren_plotting_tools as spt
import unyt
from scipy import stats


scatterplot_kwargs = {"alpha": 0.6, "s": 2, "marker": "o", "linewidth": 0.0}

# Read in cmdline arg: Are we plotting only one snapshot, or all?
# WARNING: The reference solution is comparable with snapshot_500 only
plot_all = False
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True
    snapnr = -1


plot_SPHM1RT = False
plot_Helium = False

snapshot_base = "output"


def plot_solution(filename):

    # Read in data first
    print("working on", filename)

    data = swiftsimio.load(filename)
    meta = data.metadata
    boxsize = meta.boxsize
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    xstar = data.stars.coordinates
    xpart = data.gas.coordinates
    dxp = xpart - xstar
    r = np.sqrt(np.sum(dxp ** 2, axis=1))
    r = r / (0.5 * boxsize[0])

    # Get mass fractions
    imf = spt.get_imf(scheme, data)
    xH = imf.HI + imf.HII
    xHI = imf.HI / xH
    xHII = imf.HII / xH

    # get temperature
    mu = spt.mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
    T = spt.gas_temperature(data.gas.internal_energies, mu, meta.gas_gamma)

    fig = plt.figure(figsize=(10, 5.5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.scatter(r, imf.HI, **scatterplot_kwargs, label="HI")
    ax1.scatter(r, imf.HII, **scatterplot_kwargs, label="HII")
    if plot_Helium:
        ax1.scatter(r, imf.HeI, **scatterplot_kwargs, label="HeI")
        ax1.scatter(r, imf.HeII, **scatterplot_kwargs, label="HeII")
        ax1.scatter(r, imf.HeIII, **scatterplot_kwargs, label="HeIII")
    ax1.set_ylim(5e-7, 1.5)

    ax2.scatter(r, T, **scatterplot_kwargs)

    ax1.set_xlabel("r / (L / 2)")
    ax1.set_ylabel("Hydrogen Fractions")
    ax1.set_xlim(0.0, 1.4)
    ax1.set_yscale("log")
    ax1.legend()

    ax2.set_xlabel("r / (L / 2)")
    ax2.set_ylabel("Temperature [K]")
    ax2.set_xlim(0.0, 1.0)
    ax2.set_yscale("log")

    plt.tight_layout()
    figname = filename[:-5]
    figname += ".png"

    if plot_SPHM1RT:
        figname = filename[8:-5]
        figname += "-SPHM1RT.png"

    plt.savefig(figname, dpi=200)
    plt.close()


if __name__ == "__main__":
    snaplist = spt.get_snapshot_list(snapshot_base, plot_all, snapnr)
    for f in snaplist:
        plot_solution(f)
