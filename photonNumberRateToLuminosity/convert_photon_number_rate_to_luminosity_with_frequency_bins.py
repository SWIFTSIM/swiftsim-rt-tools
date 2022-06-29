#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# Given a photon number and a minimal frequency, compute the corresponding
# luminosity in units of stellar luminosities assuming a blackbody spectrum
# for given frequency bins.
#
# Given the spectrum, we first compute the average ionizing photon energy for
# each frequency bin `i` by integrating the energy density and the number
# density from the lowest ionizing frequency to infinity (or rather until a
# very high maximal value):
#
#   <E_photon>_i = [ \int_\nu J(\nu) d\nu ] / [ \int_\nu J(\nu) / (h \nu) d\nu ]
#
# with the integration boundaries [\nu_{i, min}, \nu_{i, max}] which are
# provided by the user.
#
# Then we need to compute the fraction `f` of ionizing photons in each bin:
#
#   f_i = [ \int_{\nu_{i,min}^{\nu_{i, max}} J(\nu) / (h \nu) d\nu r ] /
#         [ \int_{\nu_{min}}^\infty J(\nu) / (h \nu) d\nu ]
#
# We finally get the luminosity L from the photon number rate Ndot as
#
#   L_i = <E_photon>_i * f_i * Ndot
#
# The final computation then simplifies to
#
#   L_i = [ \int_{\nu_{i,min}^{\nu_{i, max}} J(\nu)  d\nu r ] /
#         [ \int_{\nu_{min}}^\infty J(\nu) / (h \nu) d\nu ] * Ndot
# -----------------------------------------------------------------------------

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
integration_limit = 100 * peak_frequency
print("peak of the blackbody spectrum: {0:10.3e} [Hz]".format(peak_frequency))

# Get the full integral over the photon number density
number_density_integral_tot, nerr = integrate.quad(
    B_nu_over_h_nu, frequency_bins[0], integration_limit, args=(T, kB, h_planck, c)
)

Lsum = 0

for f in range(len(frequency_bins)):
    # Set integration ranges.
    nu_min = frequency_bins[f]
    if f < len(frequency_bins) - 1:
        nu_max = frequency_bins[f + 1]
    else:
        # A potential problem is that the integral is non-convergent
        # for the methods used. So instead of integrating to infinity,
        # peak an upper limit in units of the peak frequency.
        nu_max = integration_limit

    print("Bin {0:3d}: {1:10.3e} - {2:10.3e} [Hz] ".format(f, nu_min, nu_max), end=" ")

    energy_density_integral, eerr = integrate.quad(
        B_nu, nu_min, nu_max, args=(T, kB, h_planck, c)
    )
    L = energy_density_integral / number_density_integral_tot * Ndot / L_Sol

    print("Luminosity = {0:12.3e} [L_Sol]".format(L))

    Lsum += L

# Sanity check: Sum of individual luminosities must correspond to luminosity
# of a single frequency bin

#  energy_density_integral, eerr = integrate.quad(
#      B_nu, frequency_bins[0], integration_limit, args=(T, kB, h_planck, c)
#  )
#  number_density_integral, nerr = integrate.quad(
#      B_nu_over_h_nu, frequency_bins[0], integration_limit, args=(T, kB, h_planck, c)
#  )
#
#  L_single = energy_density_integral / number_density_integral * Ndot / L_Sol
#
#  print("Sanity check: sum over all bins / total luminosity integral = {0:12.3e}".format(Lsum / L_single))
