#!/usr/bin/env python3

# ------------------------------------------------
# Plot the cross section parametrization.
# ------------------------------------------------

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl

import cross_section_parametrization as cs_params


#  Plot parameters
params = {
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "font.size": 14,
    "font.family": "serif",
    "legend.fontsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    #  "xtick.major.width": 1.5,
    #  "ytick.major.width": 1.5,
    "axes.linewidth": 1.5,
    "text.usetex": True,
    "figure.subplot.left": 0.095,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.125,
    "figure.subplot.top": 0.875,
    "figure.subplot.wspace": 0.18,
    "figure.subplot.hspace": 0.06,
    "lines.markersize": 1,
    "lines.linewidth": 2.0,
    "mpl_toolkits.legacy_colorbar": False,
}
mpl.rcParams.update(params)



h_planck = 6.62606957e-27  # cm**2*g/s

# peak of 10^5K blackbody spectrum in Hz * 10 * h_planck
E_max = 5.879e15 * 10 * h_planck


E = np.linspace(1e-12, E_max, 10000)
cs = cs_params.PhotoionizationCrossSection()

erg_in_eV = 6.24150935e+11

cs_HI = np.zeros(E.shape)
cs_HeI = np.zeros(E.shape)
cs_HeII = np.zeros(E.shape)
for i, e in enumerate(E):
    cs_HI[i] = cs.get_cross_section(e, 0)
    cs_HeI[i] = cs.get_cross_section(e, 1)
    cs_HeII[i] = cs.get_cross_section(e, 2)


fig = plt.figure(figsize=(8,4))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

ax1.semilogy(E*erg_in_eV, cs_HI, label="HI")
ax1.semilogy(E*erg_in_eV,  cs_HeI, label="HeI")
ax1.semilogy(E*erg_in_eV, cs_HeII, label="HeII")
ax1.legend()
ax1.set_xlabel("E [eV]")
#  ax1.set_xlabel("E [erg]")
#  ax1.set_xlim(E[0], E[-1])

def tick_function(E):
    nu = E / h_planck
    return nu

ax2.set_xlim(E[0]/h_planck, E[-1]/h_planck)
ax2.set_xlabel("$\\nu$ [Hz]")


ax1.set_ylabel("$\\sigma(E)$")
#  plt.show()
plt.savefig("cross_sections.pdf")


