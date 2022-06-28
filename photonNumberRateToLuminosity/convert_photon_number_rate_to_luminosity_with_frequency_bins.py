#!/usr/bin/env python3

# ------------------------------------------------------
# Given a photon number rate, a minimal frequency, and
# frequency bins, compute the corresponding luminosity
# in units of stellar luminosities assuming a blackbody
# spectrum.
# ------------------------------------------------------

import numpy as np
from matplotlib import pyplot as plt
import unyt
import scipy.integrate as integrate

from blackbody import B_nu, B_nu_over_h_nu, nu_peak

# --------------------------------------------------------
# USER SETUP

# temperature for blackbody spectrum
T = 1e5  # K
# ionizing photon number emission rate you want
Ndot = 5e48  # s^-1
# define upper limits for frequency bins. We assume that
# the lowest bin is the first ionizing frequency.
frequency_bins = [3.288e15, 5.945e15, 13.157e15]  # Hz

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

peak_frequency = nu_peak(T, kB, h_planck)
print("peak of the blackbody spectrum: {0:10.3e} [Hz]".format(peak_frequency))

for f in range(len(frequency_bins)):
    # Set integration ranges.
    nu_min = frequency_bins[f]
    if f < len(frequency_bins) - 1:
        nu_max = frequency_bins[f + 1]
    else:
        # Don't integrate to infinity, use 100 times the peak
        # frequency of the blackbody spectrum
        nu_max = 100 * peak_frequency

    print("Bin {0:3d}: {1:10.3e} - {2:10.3e} [Hz] ".format(f, nu_min, nu_max), end=" ")

    energy_density_integral, eerr = integrate.quad(
        B_nu, nu_min, nu_max, args=(T, kB, h_planck, c)
    )
    number_density_integral, nerr = integrate.quad(
        B_nu_over_h_nu, nu_min, nu_max, args=(T, kB, h_planck, c)
    )
    L = energy_density_integral / number_density_integral * Ndot / L_Sol

    print("Luminosity = {0:12.3e} [L_Sol]".format(L))
