#!/usr/bin/env python3

#-------------------------------------------------------------------
# Plot a slice as generated by the `io_write_slice()` function call
#
# Usage:
#   python3 plotSlice.py dataFile.bin_slice_z=X.dat
#-------------------------------------------------------------------


import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import os, sys

if len(sys.argv) != 2:
    print("Error, wrong usage. Run with ./plotSlice.py dataFile.bin_slice_z=X.dat")
    quit()

infile = sys.argv[1]

if not os.path.isfile(infile):
    print("Error: file", infile, "not found")
    quit()


def set_colorbar(ax, im):
    """
    Adapt the colorbar a bit for axis object <ax> and
    imshow instance <im>
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    return




data = np.loadtxt(infile, delimiter=",")

fig = plt.figure(figsize=(6,6), dpi=200)
ax1 = fig.add_subplot(111, aspect='equal')

im1 = ax1.imshow(data, origin="lower", norm=LogNorm(),)
set_colorbar(ax1, im1)

plt.tight_layout()
fname = infile[:-3] + "png"
plt.savefig(fname)

print("Finished figure", fname)

