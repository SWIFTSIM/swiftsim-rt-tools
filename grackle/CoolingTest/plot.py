#!/usr/bin/env python3


import numpy as np
import Ptools as pt


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


fig, axs = pt.subplots(2, 2)
fig.set_size_inches(8, 6)

axs[0][0].plot(Time_Myr, Temperature, **plotkwargs)
#  axs[0].scatter(Time, Temperature)
axs[0][0].semilogy()
axs[0][0].set_xlabel(r"$\rm{Time\,\,[Myr]}$")
axs[0][0].set_ylabel(r"$\rm{Temperature\,\,[K]}$")
#  axs[0][0].set_xscale("log")

axs[0][1].plot(Time_Myr, mu)
axs[0][1].set_xlabel(r"$\rm{Time\,\,[Myr]}$")
axs[0][1].set_ylabel(r"$\rm{Mean\,\,Mol.\,\,Weight}$")
#  axs[0][1].set_xscale("log")

axs[1][0].plot(Time_Myr, nHI, ls=":", label=r"$\rm{HI}$", **plotkwargs)
axs[1][0].plot(Time_Myr, nHII, ls="-.", label=r"$\rm{HII}$", **plotkwargs)
axs[1][0].plot(Time_Myr, nHeI, ls=":", label=r"$\rm{HeI}$", **plotkwargs)
axs[1][0].plot(Time_Myr, nHeII, ls="-.", label=r"$\rm{HeII}$", **plotkwargs)
axs[1][0].plot(Time_Myr, nHeIII, ls="--", label=r"$\rm{HeIII}$", **plotkwargs)
axs[1][0].plot(Time_Myr, ne, ":", label=r"$\rm{n_e}$", **plotkwargs)
axs[1][0].plot(Time_Myr, n, label=r"$\rm{Tot}$", alpha=1)

axs[1][0].set_xlabel(r"$\rm{Time\,\,[Myr]}$")
axs[1][0].set_ylabel(r"$\rm{Number\,\,Densities}$")
#  axs[1][0].set_xscale("log")
axs[1][0].legend()

Xtot = XHI + XHII + XHeI + XHeII + XHeIII
axs[1][1].plot(Time_Myr, XHI, ls=":", label=r"$\rm{XHeI}$", **plotkwargs)
axs[1][1].plot(Time_Myr, XHII, ls="-.", label=r"$\rm{XHeII}$", **plotkwargs)
axs[1][1].plot(Time_Myr, XHeI, ls=":", label=r"$\rm{XHeI}$", **plotkwargs)
axs[1][1].plot(Time_Myr, XHeII, ls="-.", label=r"$\rm{XHeII}$", **plotkwargs)
axs[1][1].plot(Time_Myr, XHeIII, ls="--", label=r"$\rm{XHeIII}$", **plotkwargs)
axs[1][1].plot(Time_Myr, Xtot, ls="--", label="total", **plotkwargs)
axs[1][1].set_xlabel(r"$\rm{Time\,\,[Myr]}$")
axs[1][1].set_ylabel(r"$\rm{HeI,HeII,HeIII\,\,mass fraction}$")
axs[1][1].legend()
axs[1][1].grid()
#  axs[1][1].set_xscale("log")


pt.tight_layout()
pt.savefig("rt_output.png", dpi=300)
