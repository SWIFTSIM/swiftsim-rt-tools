#!/usr/bin/env python3

# ----------------------------------------------------
# Given a photon number and a minimal frequency,
# compute the corresponding luminosity in units of
# stellar luminosities assuming a blackbody spectrum
# ----------------------------------------------------

import numpy as np
from matplotlib import pyplot as plt
import unyt
import scipy.integrate as integrate

from blackbody import B_nu, B_nu_over_h_nu, nu_peak
import constants

# temperature for blackbody spectrum
T = 1e5  # K
# photon number emission rate you want
Ndot = 5e48  # s^-1
# lowest frequency
nu_min = 3.288e15  # Hz

# Use cgs values only, no unyts.
kB = constants.kB.v
h_planck = constants.h_planck.v
c = constants.c.v
L_Sol = (1 * unyt.Lsun).to("erg/s").v


nu_max = 10 * nu_peak(T, kB, h_planck)
print("nu_min: {0:12.3e} nu_peak: {1:12.3e} [Hz]".format(nu_min, nu_max))

E_scipy_numax = integrate.quad(B_nu, nu_min, nu_max, args=(T, kB, h_planck, c))
N_scipy_numax = integrate.quad(
    B_nu_over_h_nu, nu_min, nu_max, args=(T, kB, h_planck, c)
)
print(
    "Luminosity up to  10 nu_peak: {0:12.6e} [L_Sol]".format(
        E_scipy_numax[0] / N_scipy_numax[0] * Ndot / L_Sol
    )
)

nu_max = 100 * nu_peak(T, kB, h_planck)
E_scipy_numax = integrate.quad(B_nu, nu_min, nu_max, args=(T, kB, h_planck, c))
N_scipy_numax = integrate.quad(
    B_nu_over_h_nu, nu_min, nu_max, args=(T, kB, h_planck, c)
)
print(
    "Luminosity up to 100 nu_peak: {0:12.6e} [L_Sol]".format(
        E_scipy_numax[0] / N_scipy_numax[0] * Ndot / L_Sol
    )
)
