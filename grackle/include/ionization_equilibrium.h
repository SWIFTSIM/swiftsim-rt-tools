#ifndef IONIZATION_EQUILIBRIUM_H
#define IONIZATION_EQUILIBRIUM_H

#include "constants.h"
#include <grackle.h> /* needed for float type */

/*********************************************************************
/ Function to get the number densities for primordial species as a
/ function of the temperature, assuming ionization equilibrium
/ (see Katz et al. 1996, Table 2)
*********************************************************************/

gr_float A_Hp(gr_float T) {
  /* Recombination rate for H+ in units of cm^3 s^-1 */
  gr_float T3 = T / 1e3;
  gr_float T6 = T / 1e6;
  /* TODO: check this. Yves had 6.3e-11, Katz Table 2 has 8.4e-11 */
  /* return 6.30e-11 / sqrt(T) / pow(T3, 0.2) / (1 + pow(T6, 0.7)); */
  return 8.40e-11 / sqrt(T) / pow(T3, 0.2) / (1 + pow(T6, 0.7));
}

gr_float A_d(gr_float T) {
  /* dielectric recombination of He+ in units of cm^3 s^-1 */
  return 1.9e-3 / pow(T, 1.50) * exp(-470000.0 / T) *
         (1.0 + 0.30 * exp(-94000.0 / T));
}

gr_float A_Hep(gr_float T) {
  /* Recombination rate for He+ */
  /* in units of cm^3 s^-1 */
  return 1.5e-10 / pow(T, 0.6353) + A_d(T);
}

gr_float A_Hepp(gr_float T) {
  /* Recombination rate for He++ */
  /* in units of cm^3 s^-1 */
  gr_float T3 = T * 1e-3;
  gr_float T6 = T * 1e-6;
  return 3.36e-10 / sqrt(T) / pow(T3, 0.2) / (1.0 + pow(T6, 0.7));
}

gr_float G_H0(gr_float T) {
  /* (collisional) ionization rate for H0 */
  /* in units of cm^3 s^-1 */
  gr_float T5 = T / 1e5;
  return 1.17e-10 * sqrt(T) * exp(-157809.1 / T) / (1.0 + sqrt(T5));
}

gr_float G_He0(gr_float T) {
  /* in units of cm^3 s^-1 */
  gr_float T5 = T / 1e5;
  return 2.38e-11 * sqrt(T) * exp(-285335.4 / T) / (1.0 + sqrt(T5));
}

gr_float G_Hep(gr_float T) {
  /* in units of cm^3 s^-1 */
  gr_float T5 = T / 1e5;
  return 5.68e-12 * sqrt(T) * exp(-631515.0 / T) / (1.0 + sqrt(T5));
}

void ionization_equilibrium_calculate_densities(gr_float T, gr_float nH,
                                                gr_float X, gr_float *nH0,
                                                gr_float *nHp, gr_float *nHe0,
                                                gr_float *nHep, gr_float *nHepp,
                                                gr_float *ne) {
  /* Calculates densities of all species assuming
   * ionization equilibrium. This function is used to
   * compute initial conditions. */

  /* here, all returned quantities have the same units, */
  /* particles per unit volumes */

  gr_float Y = 1.0 - X;
  gr_float y = Y / (4.0 - 4.0 * Y);

  if (T < 5e3) {
    /* Below 5000K, we just go with fully neutral gas.
     * Not only is this the case, but some divisions
     * would otherwise not be safe and lead to problems
     * if we don't exception handle it here. */
    *nH0 = nH;
    *nHp = TINY_NUMBER;
    *nHe0 = y * nH;
    *nHep = TINY_NUMBER;
    *nHepp = TINY_NUMBER;
    *ne = *nHp;
  } else {

    /* Katz 1996 eq. 33 - 38 */
    /* Note: We assume all photoionization rates to be zero. */
    *nH0 = nH * A_Hp(T) / (A_Hp(T) + G_H0(T));
    *nHp = nH - *nH0;
    *nHep =
        y * nH / (1. + (A_Hep(T) + A_d(T)) / G_He0(T) + G_Hep(T) / A_Hepp(T));
    *nHe0 = *nHep * (A_Hep(T) + A_d(T)) / G_He0(T);
    *nHepp = *nHep * G_Hep(T) / A_Hepp(T);
    *ne = *nHp + *nHep + 2. * *nHepp;
  }

  /* Don't use exact zeros, the equation solver might not like it too much */
  *nH0 = fmax(*nH0, TINY_NUMBER);
  *nHp = fmax(*nHp, TINY_NUMBER);
  *nHe0 = fmax(*nHe0, TINY_NUMBER);
  *nHep = fmax(*nHep, TINY_NUMBER);
  *nHepp = fmax(*nHepp, TINY_NUMBER);
  *ne = fmax(*ne, TINY_NUMBER);
}

#endif /* IONIZATION_EQUILIBRIUM_H */
