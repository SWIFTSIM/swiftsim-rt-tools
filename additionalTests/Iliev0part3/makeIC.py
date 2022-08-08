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

from swiftsimio import Writer
from swiftsimio.units import cosmo_units

import unyt
import numpy as np
import h5py

us = cosmo_units

mass_units = us["mass"]
time_units = us["time"]
length_units = us["length"]
temperature_units = us["temperature"]
velocity_units = length_units / time_units
energy_units = mass_units * velocity_units ** 2
internal_energy_units = velocity_units ** 2

tiny_number = 1e-20

# number of particles in each dimension
n_p = 10
nparts = n_p ** 3
# filename of ICs to be generated
outputfilename = "ilievTest0part3.hdf5"
# adiabatic index
gamma = 5.0 / 3.0
# total hydrogen mass fraction
XH = 1.0
# total helium mass fraction
XHe = tiny_number
# boxsize
boxsize = (1.0 * unyt.kpc).to(length_units)
# initial gas temperature
initial_temperature = 1e2 * temperature_units
# particle mass
# take 1 amu/cm^3
pmass = (1.0 * unyt.atomic_mass_unit.to("g") / unyt.cm ** 3) * (boxsize ** 3 / nparts)
pmass = pmass.to(mass_units)

# -----------------------------------------------


def internal_energy(T, mu):
    """
    Compute the internal energy of the gas for a given
    temperature and mean molecular weight
    """
    # Using u = 1 / (gamma - 1) * p / rho
    #   and p = N/V * kT = rho / (mu * m_u) * kT

    u = unyt.boltzmann_constant * T / (gamma - 1) / (mu * unyt.atomic_mass_unit)
    return u


def mean_molecular_weight(XH0, XHp, XHe0, XHep, XHepp):
    """
    Determines the mean molecular weight for given 
    mass fractions of
        hydrogen:   XH0
        H+:         XHp
        He:         XHe0
        He+:        XHep
        He++:       XHepp

    returns:
        mu: mean molecular weight [in atomic mass units]
        NOTE: to get the actual mean mass, you still need
        to multiply it by m_u, as is tradition in the formulae
    """

    # 1/mu = sum_j X_j / A_j * (1 + E_j)
    # A_H    = 1, E_H    = 0
    # A_Hp   = 1, E_Hp   = 1
    # A_He   = 4, E_He   = 0
    # A_Hep  = 4, E_Hep  = 1
    # A_Hepp = 4, E_Hepp = 2
    one_over_mu = XH0 + 2 * XHp + 0.25 * XHe0 + 0.5 * XHep + 0.75 * XHepp

    return 1.0 / one_over_mu


# assume everything is neutral initially
mu = mean_molecular_weight(XH, 0, XHe, 0.0, 0.0)
u_part = internal_energy(initial_temperature, mu).to(internal_energy_units)
pmass = pmass.to(mass_units)


xp = unyt.unyt_array(np.zeros((nparts, 3), dtype=np.float32), length_units)
dx = boxsize.to(length_units) / n_p
ind = 0
for i in range(n_p):
    x = (i + 0.5) * dx
    for j in range(n_p):
        y = (j + 0.5) * dx
        for k in range(n_p):
            z = (k + 0.5) * dx

            xp[ind] = unyt.unyt_array((x, y, z)).to(length_units)
            ind += 1

#  w = Writer(unyt.unit_systems.cgs_unit_system, boxsize, dimension=3)
w = Writer(us, boxsize, dimension=3)

w.gas.coordinates = xp
w.gas.velocities = np.zeros(xp.shape, dtype=np.float32) * (unyt.cm / unyt.s)
w.gas.masses = np.ones(nparts, dtype=np.float32) * pmass
w.gas.internal_energy = np.ones(nparts, dtype=np.float32) * u_part

# Generate initial guess for smoothing lengths based on MIPS
w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)

# If IDs are not present, this automatically generates
w.write(outputfilename)

# Now open file back up again and add RT data.
F = h5py.File(outputfilename, "r+")
header = F["Header"]
parts = F["/PartType0"]

# Create initial ionization species mass fractions.
# Assume everything is neutral initially
HIdata = np.ones((nparts), dtype=np.float32) * XH
HIIdata = np.ones((nparts), dtype=np.float32) * tiny_number
HeIdata = np.zeros((nparts), dtype=np.float32)
HeIIdata = np.zeros((nparts), dtype=np.float32)
HeIIIdata = np.zeros((nparts), dtype=np.float32)

parts.create_dataset("MassFractionHI", data=HIdata)
parts.create_dataset("MassFractionHII", data=HIIdata)
parts.create_dataset("MassFractionHeI", data=HeIdata)
parts.create_dataset("MassFractionHeII", data=HeIIdata)
parts.create_dataset("MassFractionHeIII", data=HeIIIdata)

# close up, and we're done!
F.close()
