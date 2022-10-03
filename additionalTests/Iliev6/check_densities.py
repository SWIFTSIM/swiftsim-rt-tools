#!/usr/bin/env python3

# Check the density profile.
# Usage: check_densities <filname>
# <filename> shall be a swift output hdf5 file

import swiftsimio
import unyt
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
import sys

file = sys.argv[1]

boxsize_ref = 0.8 * unyt.kpc
data = swiftsimio.load(file)

densities = data.gas.densities
xstar = data.stars.coordinates
xpart = data.gas.coordinates

dxp = xpart - xstar
r = np.sqrt(np.sum(dxp ** 2, axis=1))
r = r / boxsize_ref

nbins = int(np.cbrt(xpart.shape[0]) + 0.5)
r_bin_edges = np.linspace(0.0, 1.7, nbins + 1)
r_bin_centers = 0.5 * (r_bin_edges[:-1] + r_bin_edges[1:])
rho_binned, _, _ = stats.binned_statistic(
    r, densities, statistic="mean", bins=r_bin_edges, range=(0.0, 1.7)
)
rho_std, _, _ = stats.binned_statistic(
    r, densities, statistic="std", bins=r_bin_edges, range=(0.0, 1.7)
)

rho_binned = rho_binned * data.gas.densities.units
rho_std = rho_std * data.gas.densities.units


# expected density profile
r_0 = 91.5 * unyt.pc
n_0 = 3.2 / unyt.cm**3
rho_0 = n_0 * unyt.proton_mass
density_expect = np.zeros(r_bin_centers.shape) * rho_0
for i,ri in enumerate(r_bin_centers):
    ri_units = ri * boxsize_ref
    if ri_units <= r_0:
        density_expect[i] = rho_0
    else:
        density_expect[i] = rho_0 * (r_0 / ri_units)**2



plt.figure(figsize=(5,5), dpi=200)
plt.errorbar(r_bin_centers, rho_binned.to("g/cm**3"), yerr=rho_std.to("g/cm**3"), label="profile")
plt.semilogy(r_bin_centers, density_expect.to("g/cm**3"), label="expected profile")

figname=file[:-5]+"-density_profile.png"
plt.savefig(figname)
