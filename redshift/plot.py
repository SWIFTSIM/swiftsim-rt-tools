#!/usr/bin/env python3

import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as cf

def fit_func(x,a,b):
    return b * x**a


# File to read from
resultfile = "out.dat"

# Read in data
data = np.loadtxt(resultfile)

step = data[:,0]
a = data[:,1]
photonEnergy = data[:,2]
volume = data[:,3]
BBtemp = data[:,4]
withShiftingBB = data[0,-1]

if withShiftingBB:
    savename = "output_redshift_shifting_BB.png"
else:
    savename = "output_redshift.png"

# Functions to convert scale factors and redshift
a2z = lambda a: 1/a - 1
z2a = lambda z: 1/(z+1)

# Calculate photon energy density
photonEnergyDensity = photonEnergy / volume


# Best fit
poptE, pcovE = cf(fit_func, a, photonEnergy)
poptErho, pcovErho = cf(fit_func, a, photonEnergyDensity)
poptBB, pcovBB = cf(fit_func, a, BBtemp)

# Create figure and axes
fig = plt.figure(dpi=300, figsize=(15,5))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

# Plot photon energy
ax1.scatter(a, photonEnergy, label="Simulation")
ax1.plot(a, fit_func(a, *poptE), c="r", label=f"Fit $\propto a^{{{poptE[0]:.3f}}}$")
secax = ax1.secondary_xaxis("top", functions=(a2z, z2a))
secax.set_xlabel("Redshift")

ax1.yaxis.get_offset_text().set_position((-0.05, 1))
ax1.set_xlabel("Scale factor")
ax1.set_ylabel("Total photon energy [IU]")
ax1.set_title("Total energy")
ax1.legend()

# Plot photon energy density
ax2.scatter(a, photonEnergyDensity, label="Simulation")
ax2.plot(a, fit_func(a, *poptErho), c="r", label=f"Fit $\propto a^{{{poptErho[0]:.3f}}}$")
secax = ax2.secondary_xaxis("top", functions=(a2z, z2a))
secax.set_xlabel("Redshift")

ax2.yaxis.get_offset_text().set_position((-0.05, 1))
ax2.set_xlabel("Scale factor")
ax2.set_ylabel("Total photon energy density [IU]")
ax2.set_title("Energy density")
ax2.legend()

# Plot black body temperature
ax3.set_xlabel("Scale factor")
ax3.set_ylabel("T [K]")
ax3.set_title("Black body temperature")
secax = ax3.secondary_xaxis("top", functions=(a2z, z2a))
secax.set_xlabel("Redshift")
ax3.yaxis.get_offset_text().set_position((-0.05, 1))

if withShiftingBB:
    ax3.plot(a, fit_func(a, *poptBB), label=f"$\propto a^{{{poptBB[0]:.3f}}}$")
    ax3.legend()
else:
    ax3.plot(a, BBtemp)
    ax3.set_ylim(0.8 * BBtemp[0], 1.2 * BBtemp[0])
    
# Clean up and save plot
plt.tight_layout()
plt.savefig(savename, dpi=300)
