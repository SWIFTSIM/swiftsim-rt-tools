#!/usr/bin/env python3

# --------------------------------------------------------------------------
# Plot the results of the heating test. Plots temperature, mean molecular
# weight, and mass fractions of Hydrogen, and internal energy.
# --------------------------------------------------------------------------

import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

plotkwargs = {"alpha": 0.6}
mh = 1.67262171e-24  # Hydrogen mass in g

fig = plt.figure(figsize=(8, 8))
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)

colors = ["C0", "C4", "C3"]

for i, sim in enumerate(["average", "min", "max"]):

    resultfile = "heating_test-density_"+sim+".dat"

    # Read in units.
    f = open(resultfile, "r")
    firstline = f.readline()
    massline = f.readline()
    lengthline = f.readline()
    velline = f.readline()
    f.close()
    units = []
    for l in [massline, lengthline, velline]:
        before, after = l.split("used:")
        val, unit = after.split("[")
        val = val.strip()
        units.append(float(val))

    mass_units = units[0]
    length_units = units[1]
    velocity_units = units[2]
    time_units = velocity_units / length_units
    density_units = mass_units / length_units ** 3
    internal_energy_units = velocity_units**2

    # Read in all other data 
    data = np.loadtxt(resultfile)

    Time = data[:, 1]
    Time_Myr = Time * 1e-6
    #  dt = data[:, 2]
    Temperature = data[:, 3]
    mu = data[:, 4]
    tot_density = data[:, 5]  # mass density
    HI_density = data[:, 6]
    HII_density = data[:, 7]
    #  HeI_density = data[:, 8]
    #  HeII_density = data[:, 9]
    #  HeIII_density = data[:, 10]
    #  e_density = data[:, 11]
    internal_energy = data[:, 11]

    XHI = HI_density / tot_density
    XHII = HII_density / tot_density
    internal_energy_cgs = internal_energy * internal_energy_units

    ax1.semilogy(Time_Myr, Temperature, label=sim+" density", c=colors[i], **plotkwargs)

    ax2.plot(Time_Myr, mu, label=sim, c=colors[i], **plotkwargs)

    ax3.semilogy(Time_Myr, XHI, linestyle="--", c=colors[i], label=sim+" density $XH_{I}$", **plotkwargs)
    ax3.semilogy(Time_Myr, XHII, linestyle=":", c=colors[i], label=sim+" density $XH_{II}$", **plotkwargs)

    ax4.semilogy(Time_Myr, internal_energy_cgs, label=sim+" density $XH_{II}$", c=colors[i], **plotkwargs)

ax1.set_ylabel(r"$\rm{Temperature\,\,[K]}$")
ax2.set_ylabel(r"$\rm{Mean\,\,Mol.\,\,Weight}$")
ax3.set_ylabel(r"Mass Fractions")
ax4.set_ylabel(r"Internal Energy [cm$^2$/s$^2$]")

for ax in fig.axes:
    ax.legend(prop={'size': 6})
    ax.grid()
    ax.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
    #  ax.set_xscale("log")

plt.tight_layout()
plt.savefig("heating_test_output.png", dpi=200)
