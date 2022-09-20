#!/usr/bin/env python3

# -----------------------------------------------------------
# Convert the ouptut of readTest5 to a hdf5 file.
# The file contains both slices and profiles.
#
# Usage: python3 convertIliev5HDF5.py
# -----------------------------------------------------------

import os
import numpy as np
import h5py

# select codes and ages you want.
codes = ["C2Ray+Capreole", "Enzo", "Flash", "HART", "Licorice", "RH1D", "RSPH", "Zeus"]
# select output times
times = ["10Myr", "30Myr", "100Myr", "200Myr", "500Myr"]
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
quantity_units = ["1", "1", "g/s**2", "K", "1", "cm**(-3)"]

ncells = 128


def make_h5file(code):
    """
    Generate a hd5f file containing all quantities in slices
    and profiles for all ages.
    """

    h5fname = code + ".hdf5"
    hfile = h5py.File(h5fname, "w")

    slices = hfile.create_group("slices")
    profiles = hfile.create_group("profiles")

    for age in times:
        slice_age = slices.create_group(age)
        profile_age = profiles.create_group(age)

        for q, quantity in enumerate(quantities):

            slicefile = code + "_" + age + "_slice_" + quantity + "_z=0.dat"
            if not os.path.exists(slicefile):
                print("Didn't find file", slicefile)
                quit(1)

            sdata = np.loadtxt(slicefile, delimiter=",")

            dss = slice_age.create_dataset(
                quantity, (ncells, ncells), dtype="f", compression="gzip"
            )
            dss[:] = sdata[:]

            proffile = code + "_" + age + "_profile_" + quantity + ".dat"
            if not os.path.exists(proffile):
                print("Didn't find file", proffile)
                quit(1)

            pdata, pstd = np.loadtxt(proffile, unpack=True, delimiter=",")

            dsp = profile_age.create_dataset(
                quantity, (ncells,), dtype="f", compression="gzip"
            )
            dsp_std = profile_age.create_dataset(
                quantity + "_std", (ncells,), dtype="f", compression="gzip"
            )
            dsp[:] = pdata[:]
            dsp_std[:] = pstd[:]

            for ds in [dss, dsp]:
                ds.attrs["description"] = quantity_descriptions[q]
                ds.attrs["unyts"] = quantity_units[q]
            dsp_std.attrs["description"] = (
                quantity_descriptions[q] + " standard deviation"
            )
            dsp_std.attrs["unyts"] = quantity_units[q]

    # add README

    README = """
    This file contains the initial conditions for the Iliev et al. 2009 "Test 5" for
    the {0:s} code. All quantities are stored in units of cell width, where the
    original test contained 128^3 cells.
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
