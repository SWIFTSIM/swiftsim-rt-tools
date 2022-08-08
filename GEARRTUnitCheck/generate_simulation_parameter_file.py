#!/usr/bin/env python3

#--------------------------------------------------------------------
# Extract required parameters from a SWIFT output or IC file and
# generate the parameter file required by the GEARRT_unit_checks
# program.
#
# usage: python3 generate_simulation_parameter_file.py <file.hdf5>
#--------------------------------------------------------------------

import yaml
import sys, os
import swiftsimio
import numpy as np


particle_mass = 0.
min_density = 0.
max_density = 0.
average_density = 0.
smoothing_length = 0.
boxsize = 0.
npart = 0
mass_units = 0.
length_units = 0.
velocity_units = 0.
current_units = 0.
temperature_units = 0.


try:
    inputfile = sys.argv[1]

    if not (os.path.exists(inputfile)):
        print("didn't find file", inputfile)
        quit()

    if not (os.path.isfile(inputfile)):
        print(inputfile, "is not a file.")

except IndexError:
    print("You need to provide this script with a file to read")
    quit()


# Read in data now
data = swiftsimio.load(inputfile)
meta = data.metadata

boxsize = meta.boxsize.max().v # keep internal units
npart = meta.n_gas

mass_units = meta.units.mass.to("g").v
length_units = meta.units.length.to("cm").v
velocity_units = (meta.units.length / meta.units.time).to("cm/s").v
current_units = meta.units.current.to("A").v
temperature_units = meta.units.temperature.to("K").v

gas = data.gas
pm = gas.masses.v
particle_mass = (np.sum(pm) / npart)
try:
    rho = gas.density.v
except AttributeError:
    try:
        rho = gas.densities.v
    except AttributeError:
        rho = 0.

average_density = np.mean(rho)
if (average_density == 0.):
    min_density = 0.
    max_density = 0.
else:
    min_density = np.min(rho)
    max_density = np.max(rho)

try:
    sml = gas.smoothing_length.v
except AttributeError:
    try:
        sml = gas.smoothing_lengths.v
    except AttributeError:
        sml = 0.

smoothing_length = np.sum(sml) / npart
min_sml = np.min(sml)
max_sml = np.max(sml)

# dump yaml file
with open("simulation_parameters.yml", "w") as file:

    # Write header

    header = "# This file was written by generate_simulation_parameter_file.py\n"
    header += "# with data extracted from " + os.path.join(os.getcwd(), inputfile)
    header += "\n\n"
    file.write(header)

    particle_mass = float(particle_mass)
    min_density = float(min_density)
    max_density = float(max_density)
    average_density = float(average_density)
    smoothing_length = float(smoothing_length)
    boxsize = float(boxsize)
    npart = int(npart)
    mass_units = float(mass_units)
    length_units = float(length_units)
    velocity_units = float(velocity_units)
    current_units = float(current_units)
    temperature_units = float(temperature_units)

    parts_dict = {
                    "ParticleMass": particle_mass,
                    "minDensity": min_density,
                    "maxDensity": max_density,
                    "averageDensity": average_density,
                    "smoothingLength": smoothing_length
            }
    units_dict = {
                    "UnitMass_in_cgs": mass_units,
                    "UnitLength_in_cgs": length_units,
                    "UnitVelocity_in_cgs": velocity_units,
                    "UnitCurrent_in_cgs": current_units,
                    "UnitTemp_in_cgs": temperature_units
                    }

    global_dict = {
                    "boxsize": boxsize, 
                    "npart": npart
                }


    file_dict = {
                    "ParticleData": parts_dict,
                    "GlobalData": global_dict,
                    "InternalUnitSystem": units_dict,
            }

    yaml.dump(file_dict, file)
