#!/usr/bin/env python3

# -----------------------------------------------------------
# Convert the ouptut of readTest4IC to a hdf5 file.
#
# Usage: python3 convertIliev4ICtoHDF5.py
# -----------------------------------------------------------

import os
import numpy as np
import h5py

# hard-code filenames you want to read in
# this should be the data output from readTest4IC
densityfile = "Iliev4DensityIC.dat"
# this should be the same file as 'sources.dat' provided online by Iliev
starfile = "sources.dat"

# file name for resulting hdf5
outputfile = "IlievTest4ICData.hdf5"

for fname in [densityfile, starfile]:
    if not os.path.isfile(fname):
        print("Didn't find file", fname)
        quit(1)

# read in redshift
with open(densityfile) as f:
    f.readline()
    redshiftline = f.readline()
    redshiftline = redshiftline.strip()
    split = redshiftline.split()
    redshift = float(split[-1])

a = 1.0 / (1.0 + redshift)

Mpc = 3.08567758e24  # cm
m_p = 1.67262311e-24  # g
M_Sun = 1.98841586e33  # g
Myr = 3600 * 24 * 365 * 1e6  # s

ncells = 128
h = 0.7
boxsize = a * 0.5 / h  # physical Mpc
dx = boxsize / ncells  # physical Mpc

#  nmin = 2.433395e-05
#  nmax = 0.1242255


def convert_n_to_M(n):
    """
    convert number density of hydrogen to cell mass
    """
    boxsize_Mpc = boxsize
    dx_Mpc = boxsize_Mpc / ncells
    V_Mpc3 = dx_Mpc ** 3
    rho_gpcm3 = n * m_p
    mass_g = (rho_gpcm3 * Mpc ** 3) * V_Mpc3
    mass_Msun = mass_g / M_Sun

    return mass_Msun


def halo_mass_from_Ndotgamma(Ndotgamma):
    """
    Get the halo mass from the read in \dot{N}_gamma
    """
    Ndotgamma_units = Ndotgamma * 1e52
    fgamma = 250  # / atom
    Omega_b = 0.043
    Omega_0 = 0.27
    ts = 3 * Myr  # Myr

    M = Ndotgamma_units * Omega_0 * m_p * ts / (fgamma * Omega_b)
    M_Msun = M / M_Sun
    return M_Msun


# Read in Data
xc, yc, zc, n = np.loadtxt(densityfile, unpack=True)
xs, ys, zs, Ndotgamma = np.loadtxt(starfile, unpack=True)

nparts = n.shape[0]
if nparts != ncells ** 3:
    print(
        "Error: Incorrect number of cells read in. Got", nparts, "expect", ncells ** 3
    )
    quit(1)
nstars = xs.shape[0]
if nstars != 16:
    print("Error: Incorrect number of stars read in. Got", nstars, "expect", 16)
    quit(1)

# Open File
hfile = h5py.File(outputfile, "w")

# Write Units
units = hfile.create_group("Units")
units.attrs["Unit current in cgs (U_I)"] = 1.0
units.attrs["Unit temperature in cgs (U_T)"] = 1.0
units.attrs["Unit time in cgs (U_t)"] = 1.0
units.attrs["Unit mass in cgs (U_M)"] = M_Sun
units.attrs["Unit length in cgs (U_L)"] = Mpc

# Write Header
header = hfile.create_group("Header")
header.attrs["Boxsize"] = [boxsize, boxsize, boxsize]
header.attrs["Dimension"] = 3
header.attrs["Redshift"] = redshift
numpart = [nparts, 0, 0, 0, nstars, 0, 0]
header.attrs["NumPartTypes"] = 7
header.attrs["NumPart_ThisFile"] = numpart
header.attrs["NumPart_Total"] = numpart


# Create and write Particle Data
parts = hfile.create_group("PartType0")

part_positions = hfile.create_dataset("PartType0/Coordinates", (nparts, 3), dtype="f")
part_mass = hfile.create_dataset("PartType0/Masses", (nparts,), dtype="f")

pos_arr = np.zeros((nparts, 3))
mass_arr = np.zeros((nparts))

for i in range(nparts):
    pos_arr[i, 0] = (xc[i] + 0.5) * dx
    pos_arr[i, 1] = (yc[i] + 0.5) * dx
    pos_arr[i, 2] = (zc[i] + 0.5) * dx
    mass_arr[i] = convert_n_to_M(n[i])

part_positions[:] = pos_arr[:]
part_mass[:] = mass_arr[:]

# Create and write Star Data
stars = hfile.create_group("PartType4")
star_positions = hfile.create_dataset("PartType4/Coordinates", (nstars, 3), dtype="f")
star_mass = hfile.create_dataset("PartType4/Masses", (nstars,), dtype="f")

star_pos_arr = np.zeros((nstars, 3))
star_mass_arr = np.zeros((nstars))

for i in range(nstars):
    star_pos_arr[i, 0] = xs[i] * dx
    star_pos_arr[i, 1] = ys[i] * dx
    star_pos_arr[i, 2] = zs[i] * dx
    star_mass_arr[i] = halo_mass_from_Ndotgamma(Ndotgamma[i])

star_positions[:] = star_pos_arr[:]
star_mass[:] = star_mass_arr[:]


# add README

README = """
This file contains the initial conditions for the Iliev et al. 2006 "Test 4".
While it mimics a regular SWIFT/gadgetoid IC file in some parts, this is not
a valid IC file. Instead, it only stores some required data in a portable 
format. In particular, the gas and star coordinates and masses.

Some Notes:

- all values are stored in !!physical!! units, not comoving units that are shown
  in the Iliev et al. 2006 paper and others. This includes in particular the
  coordinates and boxsizes.
- gas (PartType0) masses and positions are prepared to be used directly as-is
- star (PartType4) particles' masses are not true star masses, but contain the
  halo mass that is used to determine their emission rates (see eq. 13 in Iliev
  et al. 2006)
"""

meta = hfile.create_group("Metadata")

readme = hfile.create_dataset("Metadata/README", (1,), dtype=h5py.string_dtype())
readme[0] = README


hfile.close()
