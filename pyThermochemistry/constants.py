#!/usr/bin/env python3

# -------------------------
# define physical constants
# -------------------------

import unyt

# Boltzmann constant
kB = unyt.boltzmann_constant_cgs
# atomic mass unit
m_u = unyt.atomic_mass_unit_cgs
# electron mass
m_e = unyt.electron_mass_cgs
# Planck constant
h_planck = unyt.planck_constant_cgs
# speed of light
c = unyt.speed_of_light_cgs

# Threshold at which we assume gas fully ionized
T_thresh = 1e5 * unyt.K
# adiabatic index
gamma = 5.0 / 3.0
