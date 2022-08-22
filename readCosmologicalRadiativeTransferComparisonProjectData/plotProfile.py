#!/usr/bin/env python3

# -------------------------------------------------------------------
# Plot a slice as generated by the `io_write_slice()` function call
#
# Usage:
#   python3 plotSlice.py dataFile.bin_slice_z=X.dat
# -------------------------------------------------------------------


import numpy as np
from matplotlib import pyplot as plt
import os, sys

if len(sys.argv) != 2:
    print("Error, wrong usage. Run with ./plotSlice.py dataFile.bin_profile_X.dat")
    quit()

infile = sys.argv[1]

if not os.path.isfile(infile):
    print("Error: file", infile, "not found")
    quit()


data, std = np.loadtxt(infile, delimiter=",", unpack=True)

fig = plt.figure(figsize=(8, 6), dpi=200)
ax1 = fig.add_subplot(111)

#  ax1.plot(range(1, data.shape[0]+1), data)
ax1.errorbar(range(1, data.shape[0]+1), data, yerr=std, capsize=2)

plt.tight_layout()
fname = infile[:-3] + "png"
plt.savefig(fname)

print("Finished figure", fname)
