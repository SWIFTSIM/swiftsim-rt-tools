#!/usr/bin/env python3

# --------------------------------------------------------------------------
# Plot the results of the cooling test. Plots temperature, mean molecular
# weight, number densities, and mass fractions of ionizing species.
# --------------------------------------------------------------------------

import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

# file to read from
resultfile = "out.dat"

plotkwargs = {"alpha": 0.4}


mh = 1.67262171e-24  # Hydrogen mass in g

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
dt = data[:, 2]
Temperature = data[:, 3]
mu = data[:, 4]
tot_density = data[:, 5]  # mass density
HI_density = data[:, 6]
HII_density = data[:, 7]
HeI_density = data[:, 8]
HeII_density = data[:, 9]
HeIII_density = data[:, 10]
e_density = data[:, 11]  # number density


# compute number density for all species
nHI = HI_density * density_units / mh
# in part per cc
nHII = HII_density * density_units / mh
# in part per cc
nHeI = HeI_density * density_units / (4 * mh)
# in part per cc
nHeII = HeII_density * density_units / (4 * mh)
# in part per cc
nHeIII = HeIII_density * density_units / (4 * mh)
# in part per cc. Grackle Convention: electron density is 
# electron number density multiplied by proton mass
ne = e_density * density_units / mh
# in part per cc
n = nHI + nHII + nHeI + nHeII + nHeIII

XHI = HI_density / tot_density
XHII = HII_density / tot_density
XHeI = HeI_density / tot_density
XHeII = HeII_density / tot_density
XHeIII = HeIII_density / tot_density


fig = plt.figure(figsize=(8, 6))
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)

#  ax1.plot(Time_Myr, Temperature, **plotkwargs)
ax1.plot(Time, Temperature)
ax1.semilogy()
#  ax1.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
ax1.set_xlabel(r"$\rm{Time\,\,[yr]}$")
ax1.set_ylabel(r"$\rm{Temperature\,\,[K]}$")
ax1.set_xscale("log")

#  ax2.plot(Time_Myr, mu)
#  ax2.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
ax2.plot(Time, mu)
ax2.set_xlabel(r"$\rm{Time\,\,[yr]}$")
ax2.set_ylabel(r"$\rm{Mean\,\,Mol.\,\,Weight}$")
ax2.set_xscale("log")

ax3.plot(Time, nHI, ls=":", label=r"$\rm{HI}$", **plotkwargs)
ax3.plot(Time, nHII, ls="-.", label=r"$\rm{HII}$", **plotkwargs)
ax3.plot(Time, nHeI, ls=":", label=r"$\rm{HeI}$", **plotkwargs)
ax3.plot(Time, nHeII, ls="-.", label=r"$\rm{HeII}$", **plotkwargs)
ax3.plot(Time, nHeIII, ls="--", label=r"$\rm{HeIII}$", **plotkwargs)
ax3.plot(Time, ne, ":", label=r"$\rm{n_e}$", **plotkwargs)
ax3.plot(Time, n, label=r"$\rm{Tot}$", alpha=1)
ax3.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
ax3.set_ylabel(r"$\rm{Number\,\,Densities}\ [\rm{cm}^{-3}]$")
ax3.set_xscale("log")
ax3.legend()

Xtot = XHI + XHII + XHeI + XHeII + XHeIII
ax4.plot(Time, XHI, ls=":", label=r"$\rm{XHeI}$", **plotkwargs)
ax4.plot(Time, XHII, ls="-.", label=r"$\rm{XHeII}$", **plotkwargs)
ax4.plot(Time, XHeI, ls=":", label=r"$\rm{XHeI}$", **plotkwargs)
ax4.plot(Time, XHeII, ls="-.", label=r"$\rm{XHeII}$", **plotkwargs)
ax4.plot(Time, XHeIII, ls="--", label=r"$\rm{XHeIII}$", **plotkwargs)
ax4.plot(Time, Xtot, ls="--", label="total", **plotkwargs)
ax4.set_xlabel(r"$\rm{Time\,\,[yr]}$")
ax4.set_ylabel(r"$\rm{HeI,HeII,HeIII\,\,mass fraction}$")
ax4.legend()
ax4.grid()
ax4.set_xscale("log")


plt.tight_layout()
plt.savefig("ilievTest.png", dpi=300)
