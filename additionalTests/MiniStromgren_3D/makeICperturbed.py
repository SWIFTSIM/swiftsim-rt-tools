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
# Add a single star in the center of a perdurbed uniform distribution
# The gas is set up with pure hydrogen gas.
# ---------------------------------------------------------------------

import random

import numpy as np
import unyt
from swiftsimio import Writer
from swiftsimio.units import cosmo_units

import stromgren_plotting_tools as spt

random.seed(666)

if __name__ == "__main__":

    gamma = 5.0 / 3.0

    # Number of partilces in each dimension
    n_p = 40

    # Physics setup
    unitL = unyt.Mpc
    edgelen = 16.5 * 1e-3 * unitL
    edgelen = edgelen.to(unitL)
    boxsize = np.array([1.0, 1.0, 1.0]) * edgelen
    dx = edgelen / n_p

    # get gas masses
    XH = 1.0  # hydrogen mass fraction
    XHe = 0.0  # helium mass fraction
    nH = 1e-3 * unyt.cm ** (-3)
    rho_gas = nH * unyt.proton_mass
    Mtot = rho_gas * edgelen ** 3
    mpart = Mtot / n_p ** 3
    mpart = mpart.to(cosmo_units["mass"])

    # get gas internal energy for a given temperature and composition
    T = 100 * unyt.K
    XHI, XHII, XHeI, XHeII, XHeIII = spt.get_mass_fractions(T, XH, XHe)
    mu = spt.mean_molecular_weight(XHI, XHII, XHeI, XHeII, XHeIII)
    internal_energy = spt.internal_energy(T, mu, gamma)

    xp = unyt.unyt_array(np.zeros((n_p ** 3, 3), dtype=np.float64), unitL)

    # shift particles by half width so that
    # (0.5, 0.5, 0.5) can be taken by star
    shift = 0.0
    if n_p % 2 == 0:
        shift = 0.5

    ind = 0
    delta = 0.4  # max delta in units of dx
    sign = [1, -1]

    # Generate hydro coordinates
    for i in range(n_p):
        x = (i + shift) * dx
        for j in range(n_p):
            y = (j + shift) * dx
            for k in range(n_p):
                z = (k + shift) * dx

                sig = sign[random.randint(0, 1)]
                amp = random.random()
                xp[ind, 0] = x + sig * amp * delta * dx

                sig = sign[random.randint(0, 1)]
                amp = random.random()
                xp[ind, 1] = y + sig * amp * delta * dx

                sig = sign[random.randint(0, 1)]
                amp = random.random()
                xp[ind, 2] = z + sig * amp * delta * dx
                ind += 1

    xs = unyt.unyt_array(np.zeros((1, 3), dtype=np.float64), unitL)
    xs[0, 0] = 0.5 * edgelen
    xs[0, 1] = 0.5 * edgelen
    xs[0, 2] = 0.5 * edgelen

    w = Writer(unit_system=cosmo_units, box_size=boxsize, dimension=3)

    # write particle positions and smoothing lengths
    w.gas.coordinates = xp
    w.stars.coordinates = xs
    w.gas.velocities = np.zeros(xp.shape) * (unitL / unyt.Myr)
    w.stars.velocities = np.zeros(xs.shape) * (unitL / unyt.Myr)
    w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)
    w.stars.generate_smoothing_lengths(boxsize=boxsize, dimension=3)
    # gotta use particle values for 1 star
    w.stars.smoothing_length[0] = w.gas.smoothing_length[0]

    w.gas.masses = np.ones(xp.shape[0], dtype=np.float64) * mpart
    w.stars.masses = np.ones(xs.shape[0], dtype=np.float64) * mpart

    w.gas.internal_energy = np.ones(xp.shape[0], dtype=np.float64) * internal_energy

    w.write("miniStromgrenPerturbed.hdf5")
