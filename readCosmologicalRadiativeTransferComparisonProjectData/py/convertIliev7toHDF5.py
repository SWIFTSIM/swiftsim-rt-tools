#!/usr/bin/env python3

# -----------------------------------------------------------
# Convert the ouptut of readTest7 to a hdf5 file.
# The file contains both slices and profiles. It creates a
# single file for each code, containing all five output times.
#
# Usage: python3 convertIliev5HDF5.py
# -----------------------------------------------------------

import os
import numpy as np
import h5py

# select codes and ages you want.
codes = ["C2Ray+Capreole", "Coral", "Flash", "Licorice", "RSPH", "Zeus-MP"]

# select output times
times = ["1Myr", "3Myr", "10Myr", "25Myr", "50Myr"]
# which quantities to look for in files
quantities = ["xHI", "xHII", "P", "T", "mach", "n"]
# descriptions for quantities
quantity_descriptions = [
    "neutral hydrogen fraction",
    "ionized hydrogen fraction",
    "pressure",
    "temperature",
    "mach number",
    "gas number density",
]
# units (representation for unyt) of the quantities, in order
quantity_units = ["1", "1", "g/cm/s**2", "K", "1", "cm**(-3)"]

ncells = 128


def make_h5file(code):
    """
    Generate a hd5f file containing all quantities in slices
    and profiles for all ages.
    """

    h5fname = code + ".hdf5"
    hfile = h5py.File(h5fname, "w")

    slices = hfile.create_group("slices")
    slices.attrs["description"] = "Slice along the z=L/2 plane"
    profiles = hfile.create_group("profiles")
    profiles.attrs[
        "description"
    ] = "Quantity profiles along the axis of symmetry (y = z = L/2)"

    for age in times:
        slice_age = slices.create_group(age)
        profile_age = profiles.create_group(age)

        for q, quantity in enumerate(quantities):

            slicefile = code + "_" + age + "_slice_" + quantity + "_z=64.dat"
            if not os.path.exists(slicefile):
                print("Didn't find file", slicefile)
                quit(1)

            sdata = np.loadtxt(slicefile, delimiter=",")

            dss = slice_age.create_dataset(
                quantity, (ncells, ncells), dtype="f", compression="gzip"
            )
            dss[:] = sdata[:]
            dss.attrs["description"] = quantity_descriptions[q]
            dss.attrs["unyts"] = quantity_units[q]

            proffile = code + "_" + age + "_profile_" + quantity + ".dat"
            if not os.path.exists(proffile):
                print("Didn't find file", proffile)
                quit(1)

            pdata = np.loadtxt(proffile, unpack=True)

            dsp = profile_age.create_dataset(
                quantity, (ncells,), dtype="f", compression="gzip"
            )
            dsp[:] = pdata[:]

            dsp.attrs["description"] = quantity_descriptions[q]
            dsp.attrs["unyts"] = quantity_units[q]

    # add README

    README = """
    This file contains the initial conditions for the Iliev et al. 2009 "Test 7" for
    the {0:s} code. All quantities are stored in units of cell width, where the
    original test contained 128^3 cells.
    The units of each quantity is stored as a string attribute of the respective 
    quantity.
    """.format(
        code
    )

    meta = hfile.create_group("Metadata")

    readme = hfile.create_dataset("Metadata/README", (1,), dtype=h5py.string_dtype())
    readme[0] = README

    hfile.close()


if __name__ == "__main__":

    for code in codes:
        make_h5file(code)
