#!/usr/bin/env python3

# -----------------------------------------------------------------------
# Plot a profile as generated by the `io_write_profile()` function call
#
# Usage:
#   python3 plotProfile.py dataFile.bin_slice_z=X.dat
# -----------------------------------------------------------------------


import numpy as np
from matplotlib import pyplot as plt
import os, sys

if len(sys.argv) != 2:
    print("Error, wrong usage. Run with ./plotProfile.py dataFile.bin_profile_X.dat")
    quit()

infile = sys.argv[1]

if not os.path.isfile(infile):
    print("Error: file", infile, "not found")
    quit()

# differentiate the case where extracted profiles have std deviations
has_std = True
try:
    data, std = np.loadtxt(infile, delimiter=",", unpack=True)
except ValueError:
    data = np.loadtxt(infile, unpack=True)
    has_std = False

fig = plt.figure(figsize=(8, 6), dpi=200)
ax1 = fig.add_subplot(111)

if has_std:
    ax1.errorbar(range(1, data.shape[0] + 1), data, yerr=std, capsize=2)
else:
    ax1.plot(range(1, data.shape[0] + 1), data)

plt.tight_layout()
fname = infile[:-3] + "png"
plt.savefig(fname)

print("Finished figure", fname)
