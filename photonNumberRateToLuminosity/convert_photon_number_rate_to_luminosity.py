#!/usr/bin/env python3

# ----------------------------------------------------
# Given a photon number and a minimal frequency,
# compute the corresponding luminosity in units of
# stellar luminosities assuming a blackbody spectrum
# ----------------------------------------------------

import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as integrate

from blackbody import B_nu, B_nu_over_h_nu, nu_peak

# --------------------------------------------------------
# USER SETUP

# temperature for blackbody spectrum
T = 1e5  # K
# photon number emission rate you want
Ndot = 5e48  # s^-1
# lowest frequency
nu_min = 3.288e15  # Hz

# --------------------------------------------------------

# Use cgs values only, no unyts.
#  import unyt
#  kB = unyt.boltzmann_constant_cgs # erg
#  h_planck = unyt.planck_constant_cgs
#  c = unyt.speed_of_light_cgs
#  L_Sol = (1 * unyt.Lsun).to("erg/s")
kB = 1.3806488e-16  # erg/K
h_planck = 6.62606957e-27  # cm**2*g/s
c = 29979245800.0  # cm/s
L_Sol = 3.827e33  # erg/s

frequency_peak = nu_peak(T, kB, h_planck)
print("nu_min: {0:12.3e} nu_peak: {1:12.3e} [Hz]".format(nu_min, frequency_peak))

# First, do it up to 10 times the peak frequency
nu_max = 10 * frequency_peak
energy_density_integral, eerr = integrate.quad(
    B_nu, nu_min, nu_max, args=(T, kB, h_planck, c)
)
number_density_integral, nerr = integrate.quad(
    B_nu_over_h_nu, nu_min, nu_max, args=(T, kB, h_planck, c)
)
print(
    "Luminosity up to  10 nu_peak: {0:12.6e} [L_Sol]".format(
        energy_density_integral / number_density_integral * Ndot / L_Sol
    )
)

# Now, do it up to 100 times the peak frequency and see how
# much you'd be missing if you cut off at 10
nu_max = 100 * frequency_peak
energy_density_integral, eerr = integrate.quad(
    B_nu, nu_min, nu_max, args=(T, kB, h_planck, c)
)
number_density_integral, nerr = integrate.quad(
    B_nu_over_h_nu, nu_min, nu_max, args=(T, kB, h_planck, c)
)
print(
    "Luminosity up to 100 nu_peak: {0:12.6e} [L_Sol]".format(
        energy_density_integral / number_density_integral * Ndot / L_Sol
    )
)
