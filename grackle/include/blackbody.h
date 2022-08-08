#ifndef BLACKBODY_H
#define BLACKBODY_H

#include "constants.h"

/* ------------------------------------------- */
/* Functions related to the blackbody spectrum */
/* ------------------------------------------- */

/**
 * @brief Return the specific intensity of the blackbody spectrum
 *
 * @param nu frequency at which to compute specific intensity
 * @param T temperature characterizing the spectrum
 * @param kB Boltzmann constant
 * @param h_planck Planck's constant
 * @param c speed of light
 */
double blackbody_spectrum_intensity(const double nu, const double T,
                                    const double kB, const double h_planck,
                                    const double c) {

  const double hnu = h_planck * nu;
  const double kT = kB * T;
  const double nu2 = nu * nu;
  double temp;
  if (hnu / kT < 1e-6) {
    /* prevent division by zero, use Taylor approximation */
    temp = kT;
  } else if (hnu / kT > 700.) {
    /* prevent infs */
    temp = 0.;
  } else {
    temp = 1. / (exp(hnu / kT) - 1.);
  }
  return 2. * hnu * nu2 / (c * c) * temp;
}

/**
 * Return the blackbody spectrum energy density
 *
 * @param nu frequency at which to compute specific intensity
 * @param T temperature characterizing the spectrum
 * @param kB Boltzmann constant
 * @param h_planck Planck's constant
 * @param c speed of light
 */
double blackbody_spectrum_energy_density(const double nu, const double T,
                                         const double kB, const double h_planck,
                                         const double c) {
  return 4. * M_PI / c * blackbody_spectrum_intensity(nu, T, kB, h_planck, c);
}

/**
 * Return the frequency at which the blackbody spectrum at a given tempterature
 * peaks.
 *
 * @param T temperature characterizing the spectrum
 * @param kB Boltzmann constant
 * @param h_planck Planck's constant
 */
double blackbody_peak_frequency(const double T, const double kB,
                                const double h_planck) {
  return 2.82144 * kB * T / h_planck;
}

#endif /* BLACKBODY_H */
