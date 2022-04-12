#!/usr/bin/env python3

# ----------------------------------------------------------
# Various cooling and ionization rates. Values taken from
# https://ui.adsabs.harvard.edu/abs/1996ApJS..105...19K/
# ----------------------------------------------------------

import numpy as np


class rates(object):
    def __init__(self):
        return

    def A_Hp(self, T):
        """
        Recombination rate for H+ in units of cm^3 s^-1
        """

        return (
            8.40e-11
            / np.sqrt(T)
            * (T * 1e-3) ** (-0.2)
            * 1.0
            / (1.0 + (T * 1e-6) ** 0.7)
        )

    def A_d(self, T):
        """
        Dielectronic recombination rate for He+ in units of cm^3 s^-1
        """
        return (
            1.9e-3
            / T ** 1.5
            * np.exp(-470000.0 / T)
            * (1.0 + 0.3 * np.exp(-94000.0 / T))
        )

    def A_Hep(self, T):
        """
        Recombination rate for He+ in units of cm^3 s^-1
        """
        return 1.5e-10 / T ** 0.6353

    def A_Hepp(self, T):
        """
        Recombination rate for He++ in units of cm^3 s^-1
        """
        return (
            3.36e-10
            / np.sqrt(T)
            * (T * 1e-3) ** (-0.2)
            * 1.0
            / (1.0 + (T * 1e-6) ** 0.7)
        )

    def G_H0(self, T):
        """
        collisional ionization rate for H0 in units of cm^3 s^-1
        """
        return (
            5.85e-11
            * np.sqrt(T)
            * np.exp(-157809.1 / T)
            * 1.0
            / (1.0 + np.sqrt(T * 1e-5))
        )

    def G_He0(self, T):
        """
        collisional ionization rate for He0 in units of cm^3 s^-1
        """
        return (
            2.38e-11
            * np.sqrt(T)
            * np.exp(-285335.4 / T)
            * 1.0
            / (1.0 + np.sqrt(T * 1e-5))
        )

    def G_Hep(self, T):
        """
        collisional ionization rate for He+ in units of cm^3 s^-1
        
        """
        return (
            5.68e-12
            * np.sqrt(T)
            * np.exp(-631515.0 / T)
            * 1.0
            / (1.0 + np.sqrt(T * 1e-5))
        )


thermochemistry_rates = rates()
