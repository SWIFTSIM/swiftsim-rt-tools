#!/usr/bin/env python3

# ------------------------------------------------
# Plot the cross section parametrization.
# ------------------------------------------------

import numpy as np
from matplotlib import pyplot as plt

import cross_section_parametrization as cs_params

h_planck = 6.62606957e-27  # cm**2*g/s

# peak of 10^5K blackbody spectrum in Hz * 10 * h_planck
E_max = 5.879e15 * 10 * h_planck


E = np.linspace(0.0, E_max, 10000)
cs = cs_params.PhotoionizationCrossSection()


cs_HI = np.zeros(E.shape)
cs_HeI = np.zeros(E.shape)
cs_HeII = np.zeros(E.shape)
for i, e in enumerate(E):
    cs_HI[i] = cs.get_cross_section(e, 0)
    cs_HeI[i] = cs.get_cross_section(e, 1)
    cs_HeII[i] = cs.get_cross_section(e, 2)


plt.figure()
plt.semilogy(E, cs_HI, label="HI")
plt.semilogy(E, cs_HeI, label="HeI")
plt.semilogy(E, cs_HeII, label="HeII")
plt.legend()
plt.xlabel("E [erg]")
plt.ylabel("sigma parametrization")
plt.show()
