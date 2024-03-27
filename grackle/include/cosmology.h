/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *
 * Adapted from swiftsim/src/cosmology.c for use here by Mladen Ivkovic, 2024
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <math.h>

#include "error.h"

#define COSMO_TABLE_ELEMENTS 30000
#define MY_GSL_WORKSPACE_SIZE 100000

struct cosmology {
  double Omega_cdm;
  double Omega_b;
  double Omega_l;
  double Omega_k;
  double Omega_nu;
  double Omega_r;
  double w_0;
  double w_a;
  double H_0;
};

/**
 * Compute comoving density from physical density.
 **/
double cosmo_get_comoving_density(double rho_p, double a) {

  return rho_p * a * a * a;
}

/**
 * Compute physical density from comoving density.
 **/
double cosmo_get_physical_density(double rho_c, double a) {

  return rho_c / (a * a * a);
}

/**
 * Compute comoving distance from physical distance
 **/
double cosmo_get_comoving_distance(double r_p, double a) { return r_p / a; }

/**
 * Compute physical distance from comoving distance
 **/
double cosmo_get_physical_distance(double r_c, double a) { return r_c * a; }

/**
 * Compute comoving specific internal energy from physical specific internal
 *energy.
 **/
double cosmo_get_comoving_internal_energy(double u_p, double a) {

  return u_p * (a * a);
}

/**
 * Compute physical specific internal energy from comoving specific internal
 *energy.
 **/
double cosmo_get_physical_internal_energy(double u_c, double a) {

  return u_c / (a * a);
}

/**
 * Compute comoving temperature from proper temperature
 **/
double cosmo_get_comoving_temperature(double T_p, double a) {
  /* Transforms the same as internal energy. */
  return cosmo_get_comoving_internal_energy(T_p, a);
}

/**
 * Compute physical temperature from proper temperature
 **/
double cosmo_get_physical_temperature(double T_c, double a) {
  /* Transforms the same as internal energy. */
  return cosmo_get_physical_internal_energy(T_c, a);
}

/**
 * @brief Returns the interpolated value from a table.
 *
 * Uses linear interpolation.
 *
 * @param table The table of value to interpolate from (should be of length
 * COSMO_TABLE_ELEMENTS).
 * @param x The value to interpolate at.
 * @param x_min The mininum of the range of x.
 * @param x_max The maximum of the range of x.
 */
static inline double interp_table(const double table[COSMO_TABLE_ELEMENTS],
                                  const double x, const double x_min,
                                  const double x_max) {

  const double xx =
      ((x - x_min) / (x_max - x_min)) * ((double)COSMO_TABLE_ELEMENTS);

  const int i = (int)xx;
  const int ii = i < COSMO_TABLE_ELEMENTS ? i : COSMO_TABLE_ELEMENTS - 1;

  if (ii < 1)
    return table[0] * xx;
  else
    return table[ii - 1] + (table[ii] - table[ii - 1]) * (xx - ii);
}

/**
 * @brief Computes the integral of the dark-energy equation of state
 * up to a scale-factor a.
 *
 * We follow the convention of Linder & Jenkins, MNRAS, 346, 573, 2003
 * and compute \f$ \tilde{w}(a) = \int_0^a\frac{1 + w(z)}{1+z}dz \f$.
 *
 * @param a The current scale-factor.
 * @param w0 The equation of state parameter at z=0.
 * @param wa The equation of state evolution parameter.
 */
__attribute__((const)) static inline double
w_tilde(const double a, const double w0, const double wa) {
  return (a - 1.) * wa - (1. + w0 + wa) * log(a);
}

/**
 * @brief Compute E(z)$.
 *
 * @param Omega_r The radiation density parameter \f$ \Omega_r \f$.
 * @param Omega_m The matter density parameter \f$ \Omega_m \f$.
 * @param Omega_k The curvature density parameter \f$ \Omega_k \f$.
 * @param Omega_l The cosmological constant density parameter \Omega_\Lambda
 * @param w0 The equation of state parameter at z=0.
 * @param wa The equation of state evolution parameter.
 * @param a The current scale-factor.
 */
__attribute__((const)) static inline double
E(const double Omega_r, const double Omega_m, const double Omega_k,
  const double Omega_l, const double w0, const double wa, const double a) {

  const double a_inv = 1. / a;

  return sqrt(Omega_r * a_inv * a_inv * a_inv * a_inv + /* Radiation */
              Omega_m * a_inv * a_inv * a_inv +         /* Matter */
              Omega_k * a_inv * a_inv +                 /* Curvature */
              Omega_l * exp(3. * w_tilde(a, w0, wa)));  /* Lambda */
}

/**
 * @brief Computes dt for the current cosmology.
 *
 * @param a The scale-factor of interest.
 * @param param The current #cosmology.
 */
double time_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_nu = c->Omega_nu;
  const double Omega_r = c->Omega_r + Omega_nu;
  const double Omega_m = c->Omega_cdm + c->Omega_b;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_l;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double H0 = c->H_0;

  const double a_inv = 1. / a;
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);
  const double H = H0 * E_z;

  return (1. / H) * a_inv;
}

/**
 * Convert H0 from km/s/Mpc to internal units.
 */
void cosmo_convert_H0_to_internal_units(struct cosmology *cosmo,
                                        const double time_units) {

  const double km_per_Mpc = 3.240779e-20;
  const double H0_cgs = cosmo->H_0 * km_per_Mpc;
  /* from cgs -> internal units, we divide by time_units.
   * However, H_0 is in s^-1, so we multiply instead. */
  cosmo->H_0 = H0_cgs * time_units;
}

/**
 * Compute tables a(t) and t(a) for a given cosmology.
 *
 * @param a_table: where to store a(t)
 * @param t_table: where to store t(a)
 * @param cosmo: filled out cosmology struct.
 * @param a_begin: initial scale factor
 * @param a_end: final scale factor
 * @param
 */
void cosmo_get_tables(double a_table[COSMO_TABLE_ELEMENTS],
                      double t_table[COSMO_TABLE_ELEMENTS],
                      struct cosmology *cosmo, double a_begin, double a_end) {

  if (a_begin <= 0.)
    error("a_begin must be > 0. Got=%g", a_begin);
  if (a_end <= 0.)
    error("a_end must be > 0. Got=%g", a_end);
  if (a_end <= a_begin)
    error("a_begin must be > a_end Got=%g, %g", a_begin, a_end);

  const double log_a_begin = log(a_begin);
  const double log_a_end = log(a_end);

  const double delta_a = (log_a_end - log_a_begin) / COSMO_TABLE_ELEMENTS;

  for (int i = 0; i < COSMO_TABLE_ELEMENTS; i++) {
    a_table[i] = exp(log_a_begin + delta_a * (i + 1));
  }

  /* Initalise the GSL workspace */
  gsl_integration_workspace *space =
      gsl_integration_workspace_alloc(MY_GSL_WORKSPACE_SIZE);

  double result, abserr;

  /* Integrate the time \int_{a_begin}^{a_table[i]} dt */
  gsl_function F = {&time_integrand, cosmo};
  for (int i = 0; i < COSMO_TABLE_ELEMENTS; i++) {
    gsl_integration_qag(&F, a_begin, a_table[i], 0, 1.0e-10,
                        MY_GSL_WORKSPACE_SIZE, GSL_INTEG_GAUSS61, space,
                        &result, &abserr);
    t_table[i] = result;
  }
}

/**
 * Compute the cosmic time (in internal units) between expansion
 * scale factors a1 and a2.
 *
 * @param a1 scale factor at start of step
 * @param a2 scale factor at end of step
 * @param a_begin first scale factor of simulation. Needed for table limits.
 * @param a_end final scale factor of simulation. Needed for table limits.
 *
 * @return dt: physical time between a1 and a2.
 */
double cosmo_get_dt(double a1, double a2, double a_begin, double a_end,
                    double t_table[COSMO_TABLE_ELEMENTS]) {

  const double log_a1 = log(a1);
  const double log_a2 = log(a2);
  const double log_a_begin = log(a_begin);
  const double log_a_end = log(a_end);

  /* Time between a_begin and a_1 */
  const double t1 = interp_table(t_table, log_a1, log_a_begin, log_a_end);

  /* Time between a_begin and a_2 */
  const double t2 = interp_table(t_table, log_a2, log_a_begin, log_a_end);

  return t2 - t1;
}
#endif
