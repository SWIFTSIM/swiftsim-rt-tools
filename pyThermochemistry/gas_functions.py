#!/usr/bin/env python3

# ------------------------------------
# Gas related functions
# ------------------------------------

import constants
import unyt
import numpy as np
from thermochemistry_rates import thermochemistry_rates


def mean_molecular_weight(XH0, XHp, XHe0, XHep, XHepp):
    """
    Determines the mean molecular weight for given 
    mass fractions of
        hydrogen:   XH0
        H+:         XHp
        He:         XHe0
        He+:        XHep
        He++:       XHepp

    returns:
        mu: mean molecular weight [in atomic mass units]
        NOTE: to get the actual mean mass, you still need
        to multiply it by m_u, as is tradition in the formulae
    """

    # 1/mu = sum_j X_j / A_j * (1 + E_j)
    # A_H    = 1, E_H    = 0
    # A_Hp   = 1, E_Hp   = 1
    # A_He   = 4, E_He   = 0
    # A_Hep  = 4, E_Hep  = 1
    # A_Hepp = 4, E_Hepp = 2
    one_over_mu = XH0 + 2 * XHp + 0.25 * XHe0 + 0.5 * XHep + 0.75 * XHepp

    return 1.0 / one_over_mu


def gas_temperature(u, mu):
    """
    Compute the gas temperature given the specific internal 
    energy u and the mean molecular weight mu
    """

    # Using u = 1 / (gamma - 1) * p / rho
    #   and p = N/V * kT = rho / (mu * m_u) * kT

    T = u * (constants.gamma - 1) * mu * constants.m_u / constants.kB

    return T


def internal_energy(T, mu):
    """
    Compute the internal energy of the gas for a given
    temperature and mean molecular weight
    """
    # Using u = 1 / (gamma - 1) * p / rho
    #   and p = N/V * kT = rho / (mu * m_u) * kT

    u = constants.kB * T / (constants.gamma - 1) / (mu * constants.m_u)
    return u


def internal_energy_derivative(T, mu):
    """
    Compute the derivative of the internal energy of 
    the gas with respect to temperature at a given
    temperature and mean molecular weight
    """
    # Using u = 1 / (gamma - 1) * p / rho
    #   and p = N/V * kT = rho / (mu * m_u) * kT

    dudT = constants.kB / (constants.gamma - 1) / (mu * constants.m_u)
    return dudT


def get_number_densities(Temp, XH, XHe):
    """
    Compute the number densities of all species at a given
    temperature following Katz, Hernquist, and Weinberg 1996

    Temp: temperature [unyt quantity]
    XH: total mass fraction of all hydrogen species (HI + HII)
    XHe: total mass fraction of all helium species (HeI + HeII + HeIII)
    """

    # n_H = X_H * rho_gas / m_H =
    # n_He = X_He * rho_gas / m_He = (1 - X_H) / (4 X_H) * n_H
    #      =  X_He / 4(1 - X_He) * nH = y * nH

    if XH == 0:
        nH = 0.0
        nHe = 1.0
    else:
        nH = XH
        nHe = XHe / 4

    # NOTE: This is not the ionization threshold!
    if Temp < 5000 * unyt.K:
        nH0 = nH
        nHp = 0.0
        nHe0 = nHe
        nHep = 0.0
        nHepp = 0.0

    else:

        Temp.convert_to_cgs()
        T = Temp.v
        # Recombination rate for H+ in units of cm^3 s^-1
        A_Hp = thermochemistry_rates.A_Hp(T)
        # Dielectronic recombination rate for He+ in units of cm^3 s^-1
        A_d = thermochemistry_rates.A_d(T)
        # Recombination rate for He+ in units of cm^3 s^-1
        A_Hep = thermochemistry_rates.A_Hep(T)
        # Recombination rate for He++ in units of cm^3 s^-1
        A_Hepp = thermochemistry_rates.A_Hepp(T)
        # collisional ionization rate for H0 in units of cm^3 s^-1
        G_H0 = thermochemistry_rates.G_H0(T)
        # collisional ionization rate for He0 in units of cm^3 s^-1
        G_He0 = thermochemistry_rates.G_He0(T)
        # collisional ionization rate for He+ in units of cm^3 s^-1
        G_Hep = thermochemistry_rates.G_Hep(T)

        # Katz et al. 1996 eq. 33 - 38
        # Note: We assume all photoionization rates to be zero.
        # Also, here we don't care about the actual number density, i.e.
        # about the volume, since it'll cancel out later when we compute
        # the mass fractions.

        nH0 = nH * A_Hp / (A_Hp + G_H0)
        nHp = nH - nH0
        nHep = nHe / (1.0 + (A_Hep + A_d) / G_He0 + G_Hep / A_Hepp)
        nHe0 = nHep * (A_Hep + A_d) / G_He0
        nHepp = nHep * G_Hep / A_Hepp

    # electron density
    ne = nHp + nHep + 2.0 * nHepp

    return nH0, nHp, nHe0, nHep, nHepp, ne


def get_number_densities_array(Temp, XH, XHe):
    """
    Compute the number densities of all species at a given
    temperature following Katz, Hernquist, and Weinberg 1996

    Temp: temperature [unyt array]
    XH: total mass fraction of all hydrogen species (HI + HII)
    XHe: total mass fraction of all helium species (HeI + HeII + HeIII)
    """

    # n_H = X_H * rho_gas / m_H =
    # n_He = X_He * rho_gas / m_He = (1 - X_H) / (4 X_H) * n_H
    #      =  X_He / 4(1 - X_He) * nH = y * nH

    #  if XH == 0:
    #      nH = 0.0
    #      nHe = 1.0
    #  else:
    #      nH = XH
    #      nHe = XHe / 4

    nH = np.zeros(XH.shape, dtype=float)
    nHe = np.zeros(XH.shape, dtype=float)

    mask = XH == 0
    nH[mask] = 0.0
    nHe[mask] = 1.0

    inv_mask = np.logical_not(mask)
    nH[inv_mask] = XH[inv_mask]
    nHe[inv_mask] = 0.25 * XHe[inv_mask]

    # NOTE: This is not the ionization threshold!
    nH0 = np.zeros(XH.shape, dtype=float)
    nHp = np.zeros(XH.shape, dtype=float)
    nHe0 = np.zeros(XH.shape, dtype=float)
    nHep = np.zeros(XH.shape, dtype=float)
    nHepp = np.zeros(XH.shape, dtype=float)

    neutral = Temp < 5000 * unyt.K

    nH0[neutral] = nH[neutral]
    nHp[neutral] = 0.0
    nHe0[neutral] = nHe[neutral]
    nHep[neutral] = 0.0
    nHepp[neutral] = 0.0

    Temp.convert_to_cgs()
    T = Temp.v
    ionized = np.logical_not(neutral)

    # Recombination rate for H+ in units of cm^3 s^-1
    A_Hp = thermochemistry_rates.A_Hp(T)
    # Dielectronic recombination rate for He+ in units of cm^3 s^-1
    A_d = thermochemistry_rates.A_d(T)
    # Recombination rate for He+ in units of cm^3 s^-1
    A_Hep = thermochemistry_rates.A_Hep(T)
    # Recombination rate for He++ in units of cm^3 s^-1
    A_Hepp = thermochemistry_rates.A_Hepp(T)
    # collisional ionization rate for H0 in units of cm^3 s^-1
    G_H0 = thermochemistry_rates.G_H0(T)
    # collisional ionization rate for He0 in units of cm^3 s^-1
    G_He0 = thermochemistry_rates.G_He0(T)
    # collisional ionization rate for He+ in units of cm^3 s^-1
    G_Hep = thermochemistry_rates.G_Hep(T)

    # Katz et al. 1996 eq. 33 - 38
    # Note: We assume all photoionization rates to be zero.
    # Also, here we don't care about the actual number density, i.e.
    # about the volume, since it'll cancel out later when we compute
    # the mass fractions.

    nH0[ionized] = nH[ionized] * A_Hp[ionized] / (A_Hp[ionized] + G_H0[ionized])
    nHp[ionized] = nH[ionized] - nH0[ionized]
    nHep[ionized] = nHe[ionized] / (
        1.0
        + (A_Hep[ionized] + A_d[ionized]) / G_He0[ionized]
        + G_Hep[ionized] / A_Hepp[ionized]
    )
    nHe0[ionized] = nHep[ionized] * (A_Hep[ionized] + A_d[ionized]) / G_He0[ionized]
    nHepp[ionized] = nHep[ionized] * G_Hep[ionized] / A_Hepp[ionized]

    # electron density
    ne = nHp + nHep + 2.0 * nHepp

    return nH0, nHp, nHe0, nHep, nHepp, ne


def get_mass_fractions(T, XH, XHe):
    """
    Compute the mass fractions of all species at a
    given temperature

    T: temperature
    XH: total mass fraction of all hydrogen species (HI + HII)
    XHe: total mass fraction of all helium species (HeI + HeII + HeIII)
    """

    # first get number densities
    if isinstance(XH, np.ndarray):
        nH0, nHp, nHe0, nHep, nHepp, ne = get_number_densities_array(T, XH, XHe)
    else:
        nH0, nHp, nHe0, nHep, nHepp, ne = get_number_densities(T, XH, XHe)

    # now get mass denities in units of atomic mass units
    mH0 = nH0
    mHp = nHp
    mHe0 = 4.0 * nHe0
    mHep = 4.0 * nHep
    mHepp = 4.0 * nHepp
    # neglect electron mass contributions
    #  me = ne * constants.m_e / constants.m_u

    mtot = mH0 + mHp + mHe0 + mHep + mHepp  # + me

    XH0 = mH0 / mtot
    XHp = mHp / mtot
    XHe0 = mHe0 / mtot
    XHep = mHep / mtot
    XHepp = mHepp / mtot
    #  Xe = me / mtot

    return XH0, XHp, XHe0, XHep, XHepp
