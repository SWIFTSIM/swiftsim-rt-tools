#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# Compute energy- and number density weighted cross sections assuming a
# blackbody spectrum.
# In what follows, the index `i` stands for a frequency bin, while `j`
# stands for an ionizing species (H, He, He+)
#
#  let  nl = \nu_{i, min},
#       nu = \nu_{i,max}
#
#   <sigma>_ij^E = [ \int_{nl}^{nu} J(\nu) sigma_j(\nu) d\nu ] /
#                  [ \int_{nl}^{nu} J(\nu)              d\nu ] /
#
#   <sigma>_ij^N = [ \int_{nl}^{nu} J(\nu) sigma_j(\nu) / (h \nu) d\nu ] /
#                  [ \int_{nl}^{nu} J(\nu)              / (h \nu) d\nu ] /
#
# -----------------------------------------------------------------------------

import cross_section_parametrization as cs_params
import scipy.integrate as integrate
from blackbody import B_nu, B_nu_over_h_nu, nu_peak

# -----------------------------------------------------------------------------
# USER SETUP

# temperature for blackbody spectrum
T = 1e5  # K
# define upper limits for frequency bins. We assume that
# the lowest bin is the first ionizing frequency.
frequency_bins = [0.0, 3.288e15, 5.945e15, 13.157e15]  # Hz
# NOTE: If you're only using 1 photon group, start at the lowest
# ionizing frequency, not zero!! E.g:
#  frequency_bins = [3.288e15]
# number of photon groups we are using
ngroups = len(frequency_bins)

# -----------------------------------------------------------------------------

# Use cgs values only, no unyts.
#  import unyt
#  kB = unyt.boltzmann_constant_cgs # erg
#  h_planck = unyt.planck_constant_cgs
#  c = unyt.speed_of_light_cgs
kB = 1.3806488e-16  # erg/K
h_planck = 6.62606957e-27  # cm**2*g/s
c = 29979245800.0  # cm/s

# Define functions to pass to the integrator


def energy_weighted_cross_section(nu, species, cross_section, T, kB, h_planck, c):
    """
    Return the integrand for the energy weighted cross-section.

    nu: frequency
    species: for which species to work with
    cross_section: instance of PhotoionizationCrossSection
    T: temperature in K
    kB: Boltzmann constant in cgs
    h_planck: Planck's constant in cgs
    c: speed of light in cgs
    """

    E = nu * h_planck
    sigma = cross_section.get_cross_section(E, species)
    return B_nu(nu, T, kB, h_planck, c) * sigma


def number_weighted_cross_section(nu, species, cross_section, T, kB, h_planck, c):
    """
    Return the integrand for the energy weighted cross section.

    nu: frequency
    species: for which species to work with
    cross_section: instance of PhotoionizationCrossSection
    T: temperature in K
    kB: Boltzmann constant in cgs
    h_planck: Planck's constant in cgs
    c: speed of light in cgs
    """

    E = nu * h_planck
    sigma = cross_section.get_cross_section(E, species)
    return B_nu_over_h_nu(nu, T, kB, h_planck, c) * sigma


# Functions for pretty printouts


def print_header(title):
    print("\n", title)
    print("{0:20s} |".format(" "), end=" ")
    for s in range(cs_params.nspecies):
        print("         species ={0:2d} |".format(s), end=" ")
    print()
    print("---------------------", end="-")
    for s in range(cs_params.nspecies):
        print("----------------------", end="-")
    print()
    return


def print_results(results_list):
    for g in range(ngroups):

        print("            group {0:2d} |".format(g), end=" ")
        for s in range(cs_params.nspecies):
            print("{0:20.6e} |".format(results_list[g][s]), end=" ")
        print()
    return


if __name__ == "__main__":

    peak_frequency = nu_peak(T, kB, h_planck)
    integration_limit = 100 * peak_frequency
    print("peak of the blackbody spectrum: {0:10.3e} [Hz]".format(peak_frequency))

    # lists to store results in
    cse = [[] for g in range(ngroups)]
    csn = [[] for g in range(ngroups)]

    # grab an instance of the cross section parametrization
    cs = cs_params.PhotoionizationCrossSection()

    for g in range(ngroups):
        # Set integration ranges.
        nu_min = frequency_bins[g]
        if g < len(frequency_bins) - 1:
            nu_max = frequency_bins[g + 1]
        else:
            # A potential problem is that the integral is non-convergent
            # for the methods used. So instead of integrating to infinity,
            # peak an upper limit in units of the peak frequency.
            nu_max = integration_limit

        # Compute the denominators, which are species independent

        E_nu_integral, eerr = integrate.quad(
            B_nu, nu_min, nu_max, args=(T, kB, h_planck, c)
        )
        N_nu_integral, nerr = integrate.quad(
            B_nu_over_h_nu, nu_min, nu_max, args=(T, kB, h_planck, c)
        )

        for s in range(cs_params.nspecies):

            energy_weighted_integral, eerr = integrate.quad(
                energy_weighted_cross_section,
                nu_min,
                nu_max,
                args=(s, cs, T, kB, h_planck, c),
            )
            number_weighted_integral, nerr = integrate.quad(
                number_weighted_cross_section,
                nu_min,
                nu_max,
                args=(s, cs, T, kB, h_planck, c),
            )

            cse_species = energy_weighted_integral / E_nu_integral
            cse[g].append(cse_species)
            csn_species = number_weighted_integral / N_nu_integral
            csn[g].append(csn_species)

    print_header("Energy Weighted Cross Sections")
    print_results(cse)
    print_header("Number Weighted Cross Sections")
    print_results(csn)
