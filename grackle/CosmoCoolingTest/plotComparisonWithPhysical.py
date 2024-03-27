#!/usr/bin/env python3

# --------------------------------------------------------------------------
# Make comparison plots between a run with cosmology enabled, and a run
# without cosmology enabled.
# You'll have to manually tweak the `with_cosmo` variable in main.c
# and compile and run the program twice to get the output data for both
# cases.
# --------------------------------------------------------------------------

import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

plotkwargs = {"alpha": 0.4}

mh = 1.67262171e-24  # Hydrogen mass in g


fig = plt.figure(figsize=(12, 6))
ax1 = fig.add_subplot(2, 3, 1)
ax2 = fig.add_subplot(2, 3, 2)
ax3 = fig.add_subplot(2, 3, 3)
ax4 = fig.add_subplot(2, 3, 4)
ax5 = fig.add_subplot(2, 3, 5)
ax6 = fig.add_subplot(2, 3, 6)


plotnr = 0
for file, label in [("out.dat", "comoving"), ("outPhysical.dat", "physical")]:

    # Read in units.
    f = open(file, "r")
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
    data = np.loadtxt(file)

    a = data[:, 1]
    z = data[:, 2]
    Time = data[:, 3]
    Time_Myr = Time * 1e-6
    dt = data[:, 4]
    Temperature = data[:, 5]
    mu = data[:, 6]
    tot_density = data[:, 7]  # mass density
    HI_density = data[:, 8]
    HII_density = data[:, 9]
    HeI_density = data[:, 10]
    HeII_density = data[:, 11]
    HeIII_density = data[:, 12]
    e_density = data[:, 13]  # number density
    internal_energy = data[:, 14]

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

    internal_energy = internal_energy * velocity_units ** 2

    if plotnr == 0:
        ls = "-"
    else:
        ls = "--"

    ax1.plot(Time_Myr, Temperature, ls=ls, label=label, **plotkwargs)
    ax2.plot(Time_Myr, internal_energy, ls=ls, label=label, **plotkwargs)
    ax3.plot(Time_Myr, mu, ls=ls, label=label)

    ax4.plot(Time_Myr, nHI, ls=ls, label=r"$\rm{HI}$ " + label, **plotkwargs)
    ax4.plot(Time_Myr, nHeI, ls=ls, label=r"$\rm{HeI}$ " + label, **plotkwargs)
    ax4.plot(Time_Myr, n, ls=ls, label=r"$\rm{Tot}$ " + label, alpha=1)

    ax5.plot(Time_Myr, nHII, ls=ls, label=r"$\rm{HII}$ " + label, **plotkwargs)
    ax5.plot(Time_Myr, nHeII, ls=ls, label=r"$\rm{HeII}$ " + label, **plotkwargs)
    ax5.plot(Time_Myr, nHeIII, ls=ls, label=r"$\rm{HeIII}$ " + label, **plotkwargs)
    #  ax5.plot(Time_Myr, ne, ls=ls, label=r"$\rm{n_e}$"+label, **plotkwargs)

    Xtot = XHI + XHII + XHeI + XHeII + XHeIII
    ax6.plot(Time_Myr, XHI, ls=ls, label=r"$\rm{XHI}$ " + label, **plotkwargs)
    #  ax6.plot(Time_Myr, XHII, ls=ls, label=r"$\rm{XHII}$", **plotkwargs)
    ax6.plot(Time_Myr, XHeI, ls=ls, label=r"$\rm{XHeI}$ " + label, **plotkwargs)
    #  ax6.plot(Time_Myr, XHeII, ls=ls, label=r"$\rm{XHeII}$", **plotkwargs)
    #  ax6.plot(Time_Myr, XHeIII, ls=ls, label=r"$\rm{XHeIII}$", **plotkwargs)
    ax6.plot(Time_Myr, Xtot, ls=ls, label="Tot " + label, **plotkwargs)

    plotnr += 1


ax1.semilogy()
ax1.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
ax1.set_ylabel(r"$\rm{Temperature\,\,[K]}$")
ax1.set_xscale("log")
ax1.legend()


ax2.semilogy()
ax2.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
ax2.set_ylabel(r"$\rm{Internal Energy\,}$")
ax2.set_xscale("log")
ax2.legend()


ax3.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
ax3.set_ylabel(r"$\rm{Mean\,\,Mol.\,\,Weight}$")
ax3.set_xscale("log")
ax3.legend()


ax4.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
ax4.set_ylabel(r"$\rm{Number\,\,Densities}\ [\rm{cm}^{-3}]$")
ax4.set_xscale("log")
ax4.legend()

ax5.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
ax5.set_ylabel(r"$\rm{Number\,\,Densities}\ [\rm{cm}^{-3}]$")
ax5.set_xscale("log")
ax5.legend()


ax6.set_xlabel(r"$\rm{Time\,\,[Myr]}$")
ax6.set_ylabel(r"$\rm{HeI,HeII,HeIII\,\,mass fraction}$")
ax6.legend()
ax6.grid()
ax6.set_xscale("log")

plt.tight_layout()
plt.savefig("cosmo_cooling_comparison.png", dpi=300)
