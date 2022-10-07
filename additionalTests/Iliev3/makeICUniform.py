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
# Create Iliev Test 3 ICs: Uniform gas with a uniform spherical
# dense clump in it. This script gives particles a uniform grid
# as the underlying coordinate configuration.
# ---------------------------------------------------------------------

from swiftsimio import Writer
from swiftsimio.units import cosmo_units
import stromgren_plotting_tools as spt
import unyt
import numpy as np
import h5py

gamma = 5.0 / 3.0

#  resolution = 64
resolution = 128

border_particle_width = 4

if __name__ == "__main__":

    # Set up metadata
    unitL = unyt.Mpc
    # scale box size up for boundary particles
    edgelen = 6.6 * ((resolution + 2 * border_particle_width) / resolution) * unyt.kpc
    edgelen = edgelen.to(unitL)
    boxsize = np.array([1.0, 1.0, 1.0]) * edgelen

    nparts = (resolution + 2 * border_particle_width) ** 3
    dx = 1.0 / (resolution + 2 * border_particle_width)
    xp = np.zeros((nparts, 3))
    pid = np.arange(1, nparts + 1, 1)
    h = np.ones(nparts) * 2.0 * dx
    shift = border_particle_width * dx

    ind = 0
    for i in range(resolution + 2 * border_particle_width):
        x = (i + 0.5) * dx
        for j in range(resolution + 2 * border_particle_width):
            y = (j + 0.5) * dx
            for k in range(resolution + 2 * border_particle_width):
                z = (k + 0.5) * dx

                xp[ind, 0] = x
                xp[ind, 1] = y
                xp[ind, 2] = z

                is_border = False
                is_border = (
                    is_border
                    or i < border_particle_width
                    or i >= resolution + border_particle_width
                )
                is_border = (
                    is_border
                    or j < border_particle_width
                    or j >= resolution + border_particle_width
                )
                is_border = (
                    is_border
                    or k < border_particle_width
                    or k >= resolution + border_particle_width
                )

                if is_border:
                    pid[ind] += 1000000000

                    if i < border_particle_width:
                        # mark special boundary particles
                        pid[ind] *= 2

                ind += 1

    # convert to correct units
    xp *= edgelen
    h *= edgelen

    x_clump = np.array([5.0, 3.3, 3.3]).T * unyt.kpc + shift * edgelen
    r = np.sqrt(np.sum((x_clump - xp) ** 2, axis=1))
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

    nH_clump = 0.04 * unyt.cm ** (-3)
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

    w.gas.internal_energy = internal_energy

    w.write("ilievTest3UniformGrid.hdf5")
