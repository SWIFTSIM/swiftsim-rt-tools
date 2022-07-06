#!/usr/bin/env python3

# ---------------------------------------------------------
# Parametrization of photon interaction cross sections
# for a given photon frequency nu.
#
# The parametrization is taken from Verner et al. 1996
# (ui.adsabs.harvard.edu/abs/1996ApJ...465..487V)
# via Rosdahl et al. 2013
# (ui.adsabs.harvard.edu/abs/2013MNRAS.436.2188R)
# ---------------------------------------------------------

from math import sqrt

# how many species are we dealing with?
nspecies = 3


class photoionization_cross_section:

    """
    Stores the parametrization of the photoionization cross section.
    """

    def __init__(self):
        """
        Set up the parameters.
        """
        # E_ion corresponds to [13.60, 24.59, 54.42] eV
        self.E_ion = [2.179e-11, 3.940e-11, 8.719e-11]  # erg
        # E_zero corresponds to [0.4298, 0.1361, 1.720]  eV
        self.E_zero = [6.886e-13, 2.181e-13, 2.756e-12]  # erg
        self.sigma_zero = [5.475e-14, 9.492e-16, 1.369e-14]  # cm^-2
        self.P = [2.963, 3.188, 2.963]
        self.ya = [32.88, 1.469, 32.88]
        self.yw = [0.0, 2.039, 0.0]
        self.y0 = [0.0, 0.4434, 0.0]
        self.y1 = [0.0, 2.136, 0.0]

    def get_cross_section(self, energy, species):
        """
        Returns the parametrized cross section for a given
        photon energy and ionization species.

        species 0: HI
        species 1: HeI
        species 2: HeII

        energy: expected to be in erg.
        """

        E0 = self.E_zero[species]
        E_ion = self.E_ion[species]
        y0 = self.y0[species]
        y1 = self.y1[species]
        yw = self.yw[species]
        ya = self.ya[species]
        P = self.P[species]
        sigma_0 = self.sigma_zero[species]

        if energy < E_ion:
            return 0.0

        x = energy / E0 - y0
        y = sqrt(x ** 2 + y1 ** 2)
        temp1 = y ** (0.5 * P - 5.5)
        temp2 = (1.0 + sqrt(y / ya)) ** (-P)
        cs = sigma_0 * ((x - 1.0) ** 2 + yw ** 2) * temp1 * temp2

        return cs
