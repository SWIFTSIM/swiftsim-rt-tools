#!/usr/bin/env python3


# ------------------------------------------
# Blackbody spectrum related functions.
# ------------------------------------------

import numpy as np


def B_nu(nu, T, kB, h_planck, c):
    """
    Return the blackbody energy density for
    a temperature T and frequency nu
    """
    res = 2.0 * h_planck * nu ** 3 / c ** 2 / (np.exp(h_planck * nu / kB / T) - 1.0)
    return res


def B_nu_over_h_nu(nu, T, kB, h_planck, c):
    return B_nu(nu, T, kB, h_planck, c) / (h_planck * nu)


def nu_peak(T, kB, h_planck):
    """
    Return the (approximate) frequency where the peak of the
    blackbody energy density spectrum should be
    """
    return 2.82144 * kB * T / h_planck
