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
# Create ICs for Iliev Test 4.
# ---------------------------------------------------------------------

from swiftsimio import Writer
from swiftsimio.units import cosmo_units
import stromgren_plotting_tools as spt
import unyt
import numpy as np
import h5py

gamma = 5.0 / 3.0

if __name__ == "__main__":

    ic = h5py.File("IlievTest4ICData.hdf5", "r")
    parts = ic["PartType0"]
    xp = parts["Coordinates"][:]
    mp = parts["Masses"][:]
    stars = ic["PartType4"]
    xs = stars["Coordinates"][:]
    ms = stars["Masses"][:]

    units = ic["Units"]
    unitL_cgs = units.attrs["Unit length in cgs (U_L)"]
    unitM_cgs = units.attrs["Unit mass in cgs (U_M)"]

    header = ic["Header"]
    boxsize = header.attrs["Boxsize"]

    ic.close()

    nparts = int(np.cbrt(xp.shape[0]) + 0.5)

    # Add border particles
    border_particle_width = 4
    dx = boxsize[0] / nparts
    pid = np.arange(1, mp.shape[0], 1)

    # shift all particles and stars by the border width
    shift = border_particle_width * dx
    xp += shift
    xs += shift
    boxsize += 2 * shift

    npart_border = (nparts + 2 * border_particle_width) ** 3 - nparts ** 3
    xp_border = np.zeros((npart_border, 3))
    pid_border = np.arange(1000000001, 1000000001 + npart_border + 1, 1)
    ind = 0

    for i in range(nparts + 2 * border_particle_width):
        x = (i + 0.5) * dx
        for j in range(nparts + 2 * border_particle_width):
            y = (j + 0.5) * dx
            for k in range(nparts + 2 * border_particle_width):
                z = (k + 0.5) * dx

                is_border = False
                is_border = is_border or x < shift or x > boxsize[0] - shift
                is_border = is_border or y < shift or y > boxsize[1] - shift
                is_border = is_border or z < shift or z > boxsize[2] - shift

                if is_border:
                    xp_border[ind, 0] = x
                    xp_border[ind, 1] = y
                    xp_border[ind, 2] = z
                    ind += 1

    if ind != npart_border:
        print("oh no")
        quit()

    mmean = np.mean(mp)
    mp_border = np.ones(npart_border) * mmean

    # concatenate arrays
    xp = np.concatenate((xp, xp_border), axis=0)
    mp = np.concatenate((mp, mp_border), axis=0)
    pid = np.concatenate((pid, pid_border), axis=0)

    # add unyts to everything
    unitL_units = unitL_cgs * unyt.cm
    unitL = unitL_units.to("Mpc")
    unitM_units = unitM_cgs * unyt.g
    unitM = unitM_units.to("1e10 * M_Sun")

    boxsize = np.array(boxsize) * unitL
    xp = xp * unitL
    xs = xs * unitL
    mp = mp * unitM
    ms = ms * unitM

    w = Writer(unit_system=cosmo_units, box_size=boxsize, dimension=3)

    # write particle positions and smoothing lengths
    w.gas.coordinates = xp
    w.gas.particle_ids = pid

    w.stars.coordinates = xs
    # you  got to give all particle types an ID, otherwise swiftsimio
    # will generate the IDs itself
    w.stars.particle_ids = np.arange(nparts ** 3 + 1, nparts ** 3 + 1 + xs.shape[0], 1)
    w.gas.velocities = np.zeros(xp.shape) * (unitL / unyt.Myr)
    w.stars.velocities = np.zeros(xs.shape) * (unitL / unyt.Myr)
    w.gas.masses = mp
    w.stars.masses = ms
    w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)
    hmean = w.gas.smoothing_length.mean()
    w.stars.smoothing_length = np.ones(xs.shape[0]) * hmean

    # get gas internal energy for a given temperature and composition
    T = 100 * unyt.K
    XH = 1.0  # hydrogen mass fraction
    XHe = 0.0  # helium mass fraction
    XHI, XHII, XHeI, XHeII, XHeIII = spt.get_mass_fractions(T, XH, XHe)
    mu = spt.mean_molecular_weight(XHI, XHII, XHeI, XHeII, XHeIII)
    internal_energy = spt.internal_energy(T, mu, gamma)

    w.gas.internal_energy = np.ones(xp.shape[0], dtype=np.float64) * internal_energy

    w.write("ilievTest4.hdf5")
