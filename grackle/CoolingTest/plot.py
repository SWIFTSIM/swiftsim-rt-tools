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

print("Reminder: Set units manually in plot.py if you change them in main.c")
length_units = 3.08567758e21
mass_units = 1.99848e43
density_units = mass_units / length_units ** 3
time_units = 3.08567758e16


data = np.loadtxt("out.dat")

Time = data[:, 0]
Time_Myr = Time * 1e-6
dt = data[:, 1]
Temperature = data[:, 2]
mu = data[:, 3]

tot_density = data[:, 4]  # mass density

HI_density = data[:, 5]
HII_density = data[:, 6]

HeI_density = data[:, 7]
HeII_density = data[:, 8]
HeIII_density = data[:, 9]

ne = data[:, 10]  # number density


# compute number density for all species

nHI = HI_density * density_units / (mh / mass_units)
# in part per cc
nHII = HII_density * density_units / (mh / mass_units)
# in part per cc
nHeI = HeI_density * density_units / (4 * (mh / mass_units))
# in part per cc
nHeII = HeII_density * density_units / (4 * (mh / mass_units))
# in part per cc
nHeIII = HeIII_density * density_units / (4 * (mh / mass_units))
# in part per cc
ne = ne * density_units / (mh / mass_units)
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

ax1.plot(Time_Myr, Temperature, **plotkwargs)
#  ax1.scatter(Time, Temperature)
ax1.semilogy()
ax1.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
ax1.set_ylabel(r"$\rm{Temperature\,\,[K]}$")
#  ax1.set_xscale("log")

ax2.plot(Time_Myr, mu)
ax2.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
ax2.set_ylabel(r"$\rm{Mean\,\,Mol.\,\,Weight}$")
#  ax2.set_xscale("log")

ax3.plot(Time_Myr, nHI, ls=":", label=r"$\rm{HI}$", **plotkwargs)
ax3.plot(Time_Myr, nHII, ls="-.", label=r"$\rm{HII}$", **plotkwargs)
ax3.plot(Time_Myr, nHeI, ls=":", label=r"$\rm{HeI}$", **plotkwargs)
ax3.plot(Time_Myr, nHeII, ls="-.", label=r"$\rm{HeII}$", **plotkwargs)
ax3.plot(Time_Myr, nHeIII, ls="--", label=r"$\rm{HeIII}$", **plotkwargs)
ax3.plot(Time_Myr, ne, ":", label=r"$\rm{n_e}$", **plotkwargs)
ax3.plot(Time_Myr, n, label=r"$\rm{Tot}$", alpha=1)

ax3.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
ax3.set_ylabel(r"$\rm{Number\,\,Densities}$")
#  ax3.set_xscale("log")
ax3.legend()

Xtot = XHI + XHII + XHeI + XHeII + XHeIII
ax4.plot(Time_Myr, XHI, ls=":", label=r"$\rm{XHeI}$", **plotkwargs)
ax4.plot(Time_Myr, XHII, ls="-.", label=r"$\rm{XHeII}$", **plotkwargs)
ax4.plot(Time_Myr, XHeI, ls=":", label=r"$\rm{XHeI}$", **plotkwargs)
ax4.plot(Time_Myr, XHeII, ls="-.", label=r"$\rm{XHeII}$", **plotkwargs)
ax4.plot(Time_Myr, XHeIII, ls="--", label=r"$\rm{XHeIII}$", **plotkwargs)
ax4.plot(Time_Myr, Xtot, ls="--", label="total", **plotkwargs)
ax4.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
ax4.set_ylabel(r"$\rm{HeI,HeII,HeIII\,\,mass fraction}$")
ax4.legend()
ax4.grid()
#  ax4.set_xscale("log")


plt.tight_layout()
plt.savefig("cooling_test_output.png", dpi=300)
