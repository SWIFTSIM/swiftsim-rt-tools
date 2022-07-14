#!/usr/bin/env python3

# --------------------------------------------------------------------------
# Plot the results of the cooling test. Plots temperature, mean molecular
# weight, number densities, and mass fractions of ionizing species.
# --------------------------------------------------------------------------

import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

plotkwargs = {"alpha": 0.4}
mh = 1.67262171e-24  # Hydrogen mass in g

fig = plt.figure(figsize=(8, 3))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

for sim in ["average", "min", "max"]:

    resultfile = "cooling_test-density_"+sim+".dat"

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

    # Read in all other data 
    data = np.loadtxt(resultfile)

    Time = data[:, 1]
    Time_Myr = Time * 1e-6
    #  dt = data[:, 2]
    Temperature = data[:, 3]
    mu = data[:, 4]
    #  tot_density = data[:, 5]  # mass density
    #  HI_density = data[:, 6]
    #  HII_density = data[:, 7]
    #  HeI_density = data[:, 8]
    #  HeII_density = data[:, 9]
    #  HeIII_density = data[:, 10]
    #  e_density = data[:, 11]

    ax1.plot(Time_Myr, Temperature, label=sim+" density", **plotkwargs)
    ax1.semilogy()
    ax1.set_ylabel(r"$\rm{Temperature\,\,[K]}$")
    ax1.set_xscale("log")

    ax2.plot(Time_Myr, mu, label=sim, **plotkwargs)
    ax2.set_ylabel(r"$\rm{Mean\,\,Mol.\,\,Weight}$")
    ax2.set_xscale("log")


for ax in fig.axes:
    ax.legend()
    ax.grid()
    ax.set_xlabel(r"$\rm{Time\,\,[Myr]}$")

plt.tight_layout()
plt.savefig("cooling_test_output.png", dpi=200)
