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
# The ICs enforce a 1/r^2 gas density profile.
# ---------------------------------------------------------------------

import numpy as np
import unyt
from swiftsimio import Writer
from swiftsimio.units import cosmo_units

import stromgren_plotting_tools as spt

gamma = 5.0 / 3.0

#  resolution = 64
resolution = 128

# Select how many particles in the radius of the inner flat region
# you want. You might need to play with this parameter a bit to
# get a resolution similar to the desired one.
# this also determines the particle mass.
#  N_inner = 32 # this results in roughly 128^3 particles
N_inner = 16  # this results in roughly 64^3 particles


unitL = unyt.kpc
unitM = 1e10 * unyt.M_Sun

# Do everything without units to start with.
edgelen_units = 1.6 * unyt.kpc
edgelen_units = edgelen_units.to(unitL)
edgelen = edgelen_units.to(unitL).v

r_0 = (91.5 * unyt.pc).to(unitL).v
r_0_relative = (91.5 * unyt.pc / edgelen_units).to("1").v
n_0 = 3.2 / unyt.cm ** 3

rho_0_units = n_0 * unyt.proton_mass
rho_0 = rho_0_units.to(unitM / unitL ** 3).v


if __name__ == "__main__":

    print("creating ICs")
    boxsize = np.array([1.0, 1.0, 1.0]) * edgelen_units

    # Get the inner flat region < r_0
    dx_inner = r_0_relative / N_inner

    x_inner = np.zeros((2 * N_inner) ** 3)
    y_inner = np.zeros((2 * N_inner) ** 3)
    z_inner = np.zeros((2 * N_inner) ** 3)
    ind = 0

    for i in range(-N_inner, N_inner):
        xi = dx_inner * (i + 0.5)
        for j in range(-N_inner, N_inner):
            yi = dx_inner * (j + 0.5)
            for k in range(-N_inner, N_inner):
                zi = dx_inner * (k + 0.5)

                r = np.sqrt(xi ** 2 + yi ** 2 + zi ** 2)
                if r <= r_0_relative:
                    x_inner[ind] = 0.5 + xi
                    y_inner[ind] = 0.5 + yi
                    z_inner[ind] = 0.5 + zi
                    ind += 1

    xp_inner = np.vstack((x_inner[:ind], y_inner[:ind], z_inner[:ind])).T
    xp_inner *= edgelen

    # get particle mass
    V_inner = 4.0 / 3.0 * np.pi * r_0 ** 3
    m_inner = rho_0 * V_inner
    mpart = m_inner / ind

    # with the known particle mass, we can now determine how many particles
    # of the same mass will fit in a sphere with radius r_0 -> box diagonal / 2
    r_diagonal = np.sqrt(3) * 0.5 * edgelen
    Mtot = 4.0 * np.pi * rho_0 * r_0 ** 2 * (r_diagonal - (r_0 + dx_inner * edgelen))
    N_outer = int(Mtot / mpart + 0.5)

    # get r^-2 profile. This corresponds to a mass distribusion
    # which is uniform in radius.
    rng = np.random.default_rng(31415926)
    radius = rng.uniform(low=r_0 + dx_inner * edgelen, high=r_diagonal, size=N_outer)

    # draw random numbers between 0 and 1
    cos_theta = rng.uniform(low=-1.0, high=1.0, size=N_outer)
    sin_theta = np.sqrt(1.0 - cos_theta ** 2)
    phi = 2 * np.pi * rng.uniform(size=N_outer)
    xp_outer = np.zeros((N_outer, 3))
    xp_outer[:, 0] = radius * sin_theta * np.cos(phi)
    xp_outer[:, 1] = radius * sin_theta * np.sin(phi)
    xp_outer[:, 2] = radius * cos_theta

    # move sphere to center of box
    xp_outer += 0.5 * edgelen

    mask = xp_outer[:, 0] > 0.0
    mask = np.logical_and(mask, xp_outer[:, 0] < edgelen)
    mask = np.logical_and(mask, xp_outer[:, 1] > 0.0)
    mask = np.logical_and(mask, xp_outer[:, 1] < edgelen)
    mask = np.logical_and(mask, xp_outer[:, 2] > 0.0)
    mask = np.logical_and(mask, xp_outer[:, 2] < edgelen)
    xp_outer = xp_outer[mask]

    xp = np.concatenate((xp_inner, xp_outer))
    print(
        "num parts - resolution**3",
        xp.shape[0] - resolution ** 3,
        "ratio",
        (xp.shape[0] - resolution ** 3) / resolution ** 3,
    )

    # generate particle IDs
    pid = np.arange(1000000001, 1000000001 + xp.shape[0], 1)

    # generate boundary particles
    border_particle_width = 4
    dx = 1.0 / (resolution + 2 * border_particle_width)
    # scale up edge length to fit border particles
    scale = (resolution + 2 * border_particle_width) / resolution
    edgelen_units = edgelen_units * scale
    edgelen = edgelen * scale
    boxsize = boxsize * scale

    # shift actual particles
    shift = border_particle_width * dx
    xp += shift * edgelen

    npart_border = (resolution + 2 * border_particle_width) ** 3 - resolution ** 3
    xp_border = np.zeros((npart_border, 3))
    pid_border = np.arange(1, npart_border + 1, 1)
    ind = 0

    for i in range(resolution + 2 * border_particle_width):
        x = (i + 0.5) * dx
        for j in range(resolution + 2 * border_particle_width):
            y = (j + 0.5) * dx
            for k in range(resolution + 2 * border_particle_width):
                z = (k + 0.5) * dx

                is_border = False
                is_border = is_border or x < shift or x > 1.0 - shift
                is_border = is_border or y < shift or y > 1.0 - shift
                is_border = is_border or z < shift or z > 1.0 - shift

                if is_border:
                    xp_border[ind, 0] = x
                    xp_border[ind, 1] = y
                    xp_border[ind, 2] = z
                    ind += 1

    if ind != npart_border:
        print("oh no")
        quit()

    xp_border *= edgelen

    # concatenate arraays
    xp = np.concatenate((xp, xp_border), axis=0)
    pid = np.concatenate((pid, pid_border), axis=0)

    # add units to coordinates
    xp = xp * unitL

    xs = unyt.unyt_array(
        [np.array([0.5 * edgelen, 0.5 * edgelen, 0.5 * edgelen])], unitL
    )

    # get star IDs
    sid_start = 1000000001 + xp.shape[0] + 2
    sid_stop = sid_start + 1
    sid = np.arange(sid_start, sid_stop, 1)

    # Write data down
    w = Writer(
        unit_system=cosmo_units, box_size=boxsize.to(cosmo_units["length"]), dimension=3
    )

    # write particle positions and smoothing lengths
    w.gas.coordinates = xp
    w.gas.particle_ids = pid

    w.stars.coordinates = xs
    # you  got to give all particle types an ID, otherwise swiftsimio
    # will generate the IDs itself
    w.stars.particle_ids = sid
    w.gas.velocities = np.zeros(xp.shape) * (unitL / unyt.Myr)
    w.stars.velocities = np.zeros(xs.shape) * (unitL / unyt.Myr)
    #  w.gas.smoothing_length = h
    w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)
    w.stars.smoothing_length = np.ones(xs.shape[0]) * w.gas.smoothing_length.mean()

    # get gas masses
    w.gas.masses = np.ones(xp.shape[0]) * mpart * unitM
    w.stars.masses = np.ones(xs.shape[0], dtype=np.float64) * mpart * unitM

    # get gas internal energy for a given temperature and composition
    XH = 1.0  # hydrogen mass fraction
    XHe = 0.0  # helium mass fraction
    T = 100 * unyt.K
    XHI, XHII, XHeI, XHeII, XHeIII = spt.get_mass_fractions(T, XH, XHe)
    mu = spt.mean_molecular_weight(XHI, XHII, XHeI, XHeII, XHeIII)
    internal_energy = spt.internal_energy(T, mu, gamma)

    w.gas.internal_energy = np.ones(xp.shape[0], dtype=np.float64) * internal_energy

    w.write("ilievTest6.hdf5")
