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
# dense clump in it. This script creates the ICs on a uniform grid
# with different spacings for the clump and the environment so that
# the particle masses remain constant.
# ---------------------------------------------------------------------

import numpy as np
import unyt
from swiftsimio import Writer
from swiftsimio.units import cosmo_units

import stromgren_plotting_tools as spt

gamma = 5.0 / 3.0

# select target resolution.
#  resolution = 64
resolution = 128

# how many border particles to use on each face.
border_particle_width = 4


unitL = cosmo_units["length"]
unitM = cosmo_units["mass"]

edgelen_units = 6.6 * unyt.kpc
r_clump_units = 0.8 * unyt.kpc
rho_outer_units = 2e-4 / unyt.cm ** 3 * unyt.proton_mass
rho_clump_units = 200 * rho_outer_units
x_clump_units = np.array([[5.0, 3.3, 3.3]]) * unyt.kpc
T = 8000 * unyt.K
T_clump = 40 * unyt.K

# unytless quantities
edgelen = edgelen_units.to(unitL).v
r_clump = r_clump_units.to(unitL).v
rho_outer = rho_outer_units.to(unitM / unitL ** 3).v
rho_clump = rho_clump_units.to(unitM / unitL ** 3).v
x_clump = x_clump_units.to(unitL).v


if __name__ == "__main__":

    # get particle mass based on target resolution
    Npart = resolution ** 3
    V_clump = 4.0 / 3.0 * np.pi * r_clump ** 3
    V_outer = edgelen ** 3 - V_clump
    M_clump = V_clump * rho_clump
    M_outer = V_outer * rho_outer
    Mtot = M_clump + M_outer
    mpart = Mtot / Npart

    # guess how many particles we'll use for the background
    # assume full volume is spanned with background to get number of edge particles
    Nedge = int((edgelen ** 3 * rho_outer / mpart) ** (1.0 / 3.0) + 0.5)

    xp_outer = np.zeros((Nedge ** 3, 3))
    dx_outer = edgelen / Nedge
    ind = 0

    for i in range(Nedge):
        x = (i + 0.5) * dx_outer
        for j in range(Nedge):
            y = (j + 0.5) * dx_outer
            for k in range(Nedge):
                z = (k + 0.5) * dx_outer

                r2 = (
                    (x - x_clump[0, 0]) ** 2
                    + (y - x_clump[0, 1]) ** 2
                    + (z - x_clump[0, 2]) ** 2
                )
                r = np.sqrt(r2)

                if r > r_clump:
                    xp_outer[ind, 0] = x
                    xp_outer[ind, 1] = y
                    xp_outer[ind, 2] = z
                    ind += 1

    n_outer = ind
    xp_outer = xp_outer[:n_outer]

    # guess how many particles we'll use for the background
    # assume cube of size (r_clump)^3 is spanned with equal mass particles
    Nclump = int((r_clump ** 3 * rho_clump / mpart) ** (1.0 / 3.0) + 0.5)

    xp_inner = np.zeros((8 * Nclump ** 3, 3))
    dx_inner = r_clump / Nclump
    ind = 0

    # first create a sphere around (0, 0, 0)
    for i in range(-Nclump, Nclump):
        x = (i + 0.5) * dx_inner
        for j in range(-Nclump, Nclump):
            y = (j + 0.5) * dx_inner
            for k in range(-Nclump, Nclump):
                z = (k + 0.5) * dx_inner

                r2 = x ** 2 + y ** 2 + z ** 2
                r = np.sqrt(r2)

                if r <= r_clump:
                    xp_inner[ind, 0] = x
                    xp_inner[ind, 1] = y
                    xp_inner[ind, 2] = z
                    ind += 1

    n_inner = ind
    xp_inner = xp_inner[:n_inner]

    # now shift them away from (0, 0, 0)
    xp_inner += x_clump

    print("Npart environment:       ", n_outer)
    print("Npart clump:             ", n_inner)
    print("Npart tot: (w/o boundary)", n_outer + n_inner)
    print("Npart_tot/resolution**3: ", (n_outer + n_inner) / resolution ** 3)

    xp = np.concatenate((xp_outer, xp_inner), axis=0)
    pid = np.arange(1, n_outer + n_inner + 1, 1)
    is_clump = np.ones(n_outer + n_inner, dtype=bool)
    is_clump[:n_outer] = 0

    # safety check
    DX = xp[is_clump] - x_clump
    R = np.sqrt(np.sum(DX ** 2, axis=1))
    if (R > r_clump).any():
        print("OH NOOOOOOOOOOOOO", np.count_nonzero((R > r_clump)))
        quit()

    # set up quantities for border.
    # Scale edge len to fit in border particles
    scale = (Nedge + 2 * border_particle_width) / Nedge
    dx_border = edgelen / Nedge
    edgelen *= scale
    edgelen_units = edgelen_units * scale
    # shift non-border parts
    shift = border_particle_width * dx_border
    xp += shift

    npart_border = (Nedge + 2 * border_particle_width) ** 3 - Nedge ** 3

    # get border particles now
    xp_border = np.zeros((npart_border, 3))
    pid_border = np.arange(1000000001, 1000000001 + npart_border, 1)
    is_clump_border = np.zeros(npart_border, dtype=bool)
    ind = 0

    for i in range(Nedge + 2 * border_particle_width):
        x = (i + 0.5) * dx_border
        for j in range(Nedge + 2 * border_particle_width):
            y = (j + 0.5) * dx_border
            for k in range(Nedge + 2 * border_particle_width):
                z = (k + 0.5) * dx_border

                is_border = False
                is_border = is_border or x < shift or x > edgelen - shift
                is_border = is_border or y < shift or y > edgelen - shift
                is_border = is_border or z < shift or z > edgelen - shift

                if is_border:
                    xp_border[ind, 0] = x
                    xp_border[ind, 1] = y
                    xp_border[ind, 2] = z

                    if i > border_particle_width:
                        # mark special boundary particles: add 10^7 to particles
                        # that won't "emit" radiation
                        pid_border[ind] *= 2

                    ind += 1

    if ind != npart_border:
        print("oh no", ind, npart_border)
        quit()

    # concatenate arraays
    xp = np.concatenate((xp_border, xp), axis=0)
    pid = np.concatenate((pid_border, pid), axis=0)
    is_clump = np.concatenate((is_clump_border, is_clump), axis=0)

    # convert to correct units
    xp = xp * unitL
    masses = np.ones(xp.shape[0]) * mpart * unitM

    # set up quantities assuming they're outside the clump
    XH = 1.0  # hydrogen mass fraction
    XHe = 0.0  # helium mass fraction

    # get gas internal energy for a given temperature and composition
    XHI, XHII, XHeI, XHeII, XHeIII = spt.get_mass_fractions(T, XH, XHe)
    mu = spt.mean_molecular_weight(XHI, XHII, XHeI, XHeII, XHeIII)
    u = spt.internal_energy(T, mu, gamma)

    internal_energy = np.ones(xp.shape[0]) * u

    XHI, XHII, XHeI, XHeII, XHeIII = spt.get_mass_fractions(T_clump, XH, XHe)
    mu_clump = spt.mean_molecular_weight(XHI, XHII, XHeI, XHeII, XHeIII)
    u_clump = spt.internal_energy(T_clump, mu_clump, gamma)

    internal_energy[is_clump] = u_clump.to(u.units)

    # write file
    boxsize = np.array([1.0, 1.0, 1.0]) * edgelen_units.to(unitL)

    w = Writer(unit_system=cosmo_units, box_size=boxsize, dimension=3)

    w.gas.coordinates = xp
    w.gas.particle_ids = pid
    w.gas.masses = masses
    w.gas.velocities = np.zeros(xp.shape) * (unitL / unyt.Myr)
    w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)
    w.gas.internal_energy = internal_energy

    w.write("ilievTest3EqualMasses.hdf5")
