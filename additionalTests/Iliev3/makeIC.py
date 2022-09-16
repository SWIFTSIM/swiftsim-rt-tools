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
# Add a single star in the center of a glass distribution
# The gas is set up with pure hydrogen gas.
# ---------------------------------------------------------------------

from swiftsimio import Writer
from swiftsimio.units import cosmo_units
import stromgren_plotting_tools as spt
import unyt
import numpy as np
import h5py

gamma = 5.0 / 3.0

resolution = 64
#  resolution = 128

border_particle_width = 4

if __name__ == "__main__":

    glass = h5py.File("glassCube_"+str(resolution)+".hdf5", "r")
    #  glass = h5py.File("glassCube_128.hdf5", "r")

    parts = glass["PartType0"]
    xp = parts["Coordinates"][:]
    h = parts["SmoothingLength"][:]
    glass.close()
    nparts = resolution


    # Double-check all particles for boundaries
    for i in range(3):
        mask = xp[:, i] < 0.0
        xp[mask, i] += 1.0
        mask = xp[:, i] > 1.0
        xp[mask, i] -= 1.0



    # Set up metadata
    unitL = unyt.Mpc
    # scale box size up for boundary particles
    edgelen = 6.6 * ((resolution + 2 * border_particle_width) / resolution) * unyt.kpc
    edgelen = edgelen.to(unitL)
    boxsize = np.array([1.0, 1.0, 1.0]) * edgelen

    # Add border particles
    dx = 1.0 / (nparts + 2 * border_particle_width)
    scale = nparts / (nparts + 2 * border_particle_width)
    shift = 0.5 * (1.0 - scale)
    if scale < 0:
        print("scale =", scale, "???")
        quit()

    # scale positions down
    xp *= scale
    # and shift them away from (0, 0, 0)
    xp += shift
    pid = np.arange(1, h.shape[0], 1)

    npart_border = (nparts + 2 * border_particle_width) ** 3 - nparts ** 3

    # add border particles
    xp_border = np.zeros((npart_border, 3))
    h_border = np.ones((npart_border))
    pid_border = np.arange(1000000001, 1000000001 + npart_border + 1, 1)
    ind = 0

    for i in range(nparts + 2 * border_particle_width):
        x = (i + 0.5) * dx
        for j in range(nparts + 2 * border_particle_width):
            y = (j + 0.5) * dx
            for k in range(nparts + 2 * border_particle_width):
                z = (k + 0.5) * dx

                is_border = False
                is_border = is_border or x < shift or x > 1.0 - shift
                is_border = is_border or y < shift or y > 1.0 - shift
                is_border = is_border or z < shift or z > 1.0 - shift

                if is_border:
                    xp_border[ind, 0] = x
                    xp_border[ind, 1] = y
                    xp_border[ind, 2] = z

                    if y < shift:
                        # mark special boundary particles
                        pid_border[ind] *= 2

                    ind += 1



    if ind != (nparts + 2 * border_particle_width) ** 3 - nparts ** 3:
        print("oh no")
        quit()

    h_mean = np.mean(h)
    # take estimate for sml
    h_border = h_border * h_mean

    # concatenate arraays
    xp = np.concatenate((xp, xp_border), axis=0)
    h = np.concatenate((h, h_border), axis=0)
    pid = np.concatenate((pid, pid_border), axis=0)

    # convert to correct units
    xp *= edgelen
    h *= edgelen

    x_clump = np.array([5., 3.3, 3.3]).T * unyt.kpc + shift * edgelen
    r = np.sqrt(np.sum((x_clump - xp)**2, axis=1))
    is_clump = r <= 0.8 * unyt.kpc

    # set up quantities assuming they're outside the clump
    XH = 1.0  # hydrogen mass fraction
    XHe = 0.0  # helium mass fraction
    nH = 2e-4 * unyt.cm ** (-3)

    # get particle mass assuming all volumes are equal
    rho_gas = nH * unyt.proton_mass
    Mtot = rho_gas * edgelen ** 3
    mpart = Mtot / xp.shape[0]
    mpart = mpart.to(cosmo_units["mass"])

    masses = np.ones(xp.shape[0]) * mpart

    nH_clump = 0.04 * unyt.cm**(-3)
    rho_clump = nH_clump * unyt.proton_mass
    Mtot_clump = rho_clump * edgelen ** 3
    mpart_clump = Mtot_clump / xp.shape[0]
    mpart_clump = mpart_clump.to(cosmo_units["mass"])

    masses[is_clump] = mpart_clump.to(masses.units)

    # get gas internal energy for a given temperature and composition
    T = 8000 * unyt.K
    XHI, XHII, XHeI, XHeII, XHeIII = spt.get_mass_fractions(T, XH, XHe)
    mu = spt.mean_molecular_weight(XHI, XHII, XHeI, XHeII, XHeIII)
    u = spt.internal_energy(T, mu, gamma)

    internal_energy = np.ones(xp.shape[0]) * u

    T_clump = 40 * unyt.K
    XHI, XHII, XHeI, XHeII, XHeIII = spt.get_mass_fractions(T_clump, XH, XHe)
    mu_clump = spt.mean_molecular_weight(XHI, XHII, XHeI, XHeII, XHeIII)
    u_clump = spt.internal_energy(T_clump, mu_clump, gamma)

    internal_energy[is_clump] = u_clump.to(u.units)



    # write file

    w = Writer(unit_system=cosmo_units, box_size=boxsize, dimension=3)

    w.gas.coordinates = xp
    w.gas.particle_ids = pid
    w.gas.masses = masses
    w.gas.velocities = np.zeros(xp.shape) * (unitL / unyt.Myr)
    w.gas.smoothing_length = h

    w.gas.internal_energy = np.ones(xp.shape[0], dtype=np.float64) * internal_energy

    w.write("ilievTest3.hdf5")
