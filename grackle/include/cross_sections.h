#ifndef CROSS_SECTIONS_H
#define CROSS_SECTIONS_H

/* ------------------------------------------------------ */
/* Functions related to the photoionizing cross sections. */
/* ------------------------------------------------------ */

#include "blackbody.h"
#include "constants.h"
#include <gsl/gsl_integration.h>

/* HI, HeI, HeII */
#define RT_NIONIZING_SPECIES 3
#define RT_INTEGRAL_NPOINTS 1000

/*! Struct containing the parametrized cross section parameters
 * for each photoionizing species.
 * Correct usage is to call init_photoion_cs_params_cgs(),
 * which returns a fully initialized struct. */
struct photoion_cs_parameters {
  double E_ion[RT_NIONIZING_SPECIES];
  double E_zero[RT_NIONIZING_SPECIES];
  double sigma_zero[RT_NIONIZING_SPECIES];
  double P[RT_NIONIZING_SPECIES];
  double ya[RT_NIONIZING_SPECIES];
  double yw[RT_NIONIZING_SPECIES];
  double y0[RT_NIONIZING_SPECIES];
  double y1[RT_NIONIZING_SPECIES];
};

/*! Parameters needed to compute the stellar spectra integrals -
 * the data needs to be encapsulated in a single struct to be passed
 * as an argument for GSL integrators. */
struct spectrum_integration_params {
  /* Which species are we dealing with? */
  int species;
  /* Temperature of blackbody in correct units */
  double T;
  /* Boltzmann constant in correct units */
  double kB;
  /* Planck constant in correct units */
  double h_planck;
  /* speed of light in correct units */
  double c;
  /* Values for the cross section parametrization */
  struct photoion_cs_parameters *cs_params;
};

/**
 * Initialize the parameters for the cross section computation in cgs,
 * and return a fully and correctly initialized struct.
 * The data is taken from Verner et al. 1996
 * (ui.adsabs.harvard.edu/abs/1996ApJ...465..487V) via Rosdahl et al. 2013
 * (ui.adsabs.harvard.edu/abs/2013MNRAS.436.2188R)
 */
struct photoion_cs_parameters init_photoion_cs_params_cgs(void) {

  struct photoion_cs_parameters photoion_cs_params_cgs = {
      /* E_ion =         {13.60,     24.59,     54.42}          eV */
      /* E_zero =        {0.4298,    13.61,    1.720},         eV */
      /* Please note that the value for E_0 of HeI of 0.1361 given
       * in table E1 of Rosdahl et al is wrong. It's 13.61 instead */
      /* E_ion =      */ {2.179e-11, 3.940e-11, 8.719e-11}, /* erg */
      /* E_zero =     */ {6.886e-13, 2.181e-11, 2.756e-12}, /* erg */
      /* sigma_zero = */ {5.475e-14, 9.492e-16, 1.369e-14}, /* cm^-2 */
      /* P =          */ {2.963, 3.188, 2.963},
      /* ya =         */ {32.88, 1.469, 32.88},
      /* yw =         */ {0., 2.039, 0.},
      /* y0 =         */ {0., 0.4434, 0.},
      /* y1 =         */ {0., 2.136, 0.}};

  return photoion_cs_params_cgs;
}

/**
 * Compute the parametrized cross section for a given energy and species.
 * The parametrization is taken from Verner et al. 1996
 * (ui.adsabs.harvard.edu/abs/1996ApJ...465..487V) via Rosdahl et al. 2013
 * (ui.adsabs.harvard.edu/abs/2013MNRAS.436.2188R)
 *
 * @param E energy for which to compute the cross section.
 * @param species index of species, 0 < species < RT_NIONIZING_SPECIES
 * @param params cross section parameters struct
 */
double
photoionization_cross_section(const double E, const int species,
                              const struct photoion_cs_parameters *params) {

  const double E0 = params->E_zero[species];
  const double E_ion = params->E_ion[species];
  const double y0 = params->y0[species];
  const double y1 = params->y1[species];
  const double yw = params->yw[species];
  const double ya = params->ya[species];
  const double P = params->P[species];
  const double sigma_0 = params->sigma_zero[species];

  if (E < E_ion)
    return 0.;

  const double x = E / E0 - y0;
  const double y = sqrt(x * x + y1 * y1);
  const double temp1 = pow(y, 0.5 * P - 5.5);
  const double temp2 = pow(1. + sqrt(y / ya), -P);

  return sigma_0 * ((x - 1.) * (x - 1.) + yw * yw) * temp1 * temp2;
}

/**
 * Spectrum function to be integrated.
 * This function is called by the GSL integrator.
 *
 * @param nu frequency at which to evaluate spectrum
 * @param params spectrum integration params. Needs to be of type
 * void* for GSL integrators.
 */
double spectrum_integrand(double nu, void *params) {
  struct spectrum_integration_params *pars =
      (struct spectrum_integration_params *)params;
  const double T = pars->T;
  const double kB = pars->kB;
  const double h_planck = pars->h_planck;
  const double c = pars->c;
  return blackbody_spectrum_intensity(nu, T, kB, h_planck, c);
}

/**
 * Spectrum function divided by photon energy h*nu to be integrated.
 * This function is called by the GSL integrator.
 *
 * @param nu frequency at which to evaluate spectrum
 * @param params spectrum integration params. Needs to be of type
 * void* for GSL integrators.
 */
double spectrum_over_hnu_integrand(double nu, void *params) {
  struct spectrum_integration_params *p =
      (struct spectrum_integration_params *)params;
  const double T = p->T;
  const double kB = p->kB;
  const double h_planck = p->h_planck;
  const double c = p->c;
  const double E = nu * h_planck;

  if (E > 0.) {
    const double J = blackbody_spectrum_intensity(nu, T, kB, h_planck, c);
    return J / E;
  } else {
    return 0.;
  }
}

/**
 * Spectrum times cross section function to be integrated.
 * This function is called by the GSL integrator.
 *
 * @param nu frequency at which to evaluate spectrum
 * @param params spectrum integration params. Needs to be of type
 * void* for GSL integrators.
 */
double spectrum_times_sigma_integrand(double nu, void *params) {
  struct spectrum_integration_params *p =
      (struct spectrum_integration_params *)params;
  const double h_planck = p->h_planck;
  const double E = nu * p->h_planck;
  if (E > 0.) {
    const double sigma =
        photoionization_cross_section(E, p->species, p->cs_params);
    const double T = p->T;
    const double kB = p->kB;
    const double c = p->c;
    const double J = blackbody_spectrum_intensity(nu, T, kB, h_planck, c);
    return J * sigma;
  } else {
    return 0.;
  }
}

/**
 * Spectrum times cross section divided by h*nu function to be integrated.
 * This function is called by the GSL integrator.
 *
 * @param nu frequency at which to evaluate spectrum
 * @param params spectrum integration params. Needs to be of type
 * void* for GSL integrators.
 */
double spectrum_times_sigma_over_hnu_integrand(double nu, void *params) {
  struct spectrum_integration_params *p =
      (struct spectrum_integration_params *)params;
  const double h_planck = p->h_planck;
  const double E = nu * p->h_planck;
  if (E > 0.) {
    const double sigma =
        photoionization_cross_section(E, p->species, p->cs_params);
    const double T = p->T;
    const double kB = p->kB;
    const double c = p->c;
    const double J = blackbody_spectrum_intensity(nu, T, kB, h_planck, c);
    return J * sigma / E;
  } else {
    return 0.;
  }
}

/**
 * Integrate a function from nu_start to nu_stop with GSL routines
 *
 * @param function function to integrate
 * @param nu_start lower boundary of the integral
 * @param nu_stop upper boundary of the integral
 * @param params spectrum integration params.
 * */
double
cross_sections_integrate_gsl(double (*function)(double, void *),
                             double nu_start, double nu_stop, int npoints,
                             struct spectrum_integration_params *params) {

  gsl_function F;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(npoints);
  double result, error;

  F.function = function;
  F.params = (void *)params;
  gsl_integration_qags(&F, nu_start, nu_stop, /*espabs=*/0., /*epsrel=*/1e-7,
                       npoints, w, &result, &error);

  gsl_integration_workspace_free(w);

  return result;
}

/**
 * @brief allocate and compute the averaged cross sections
 * for each photon group and ionizing species.
 *
 * Assume everything is in CGS units.
 * NOTE: all return parameters will be malloc'd in here.
 *
 * @param T_blackbody: temperature of the blackbody spectrum in K
 * @param frequency_bins: lower frequency edges for the photon group bins, in Hz
 * @return cse: energy weighted averaged cross sections.
 * @return csn: number weighted averaged cross sections.
 * @return mean_energy: mean photon energy in each photon frequency bin.
 **/
void get_cross_sections(double T_blackbody, double frequency_bins[RT_NGROUPS],
                        double **cse, double **csn,
                        double mean_energy[RT_NGROUPS]) {

  double integral_E[RT_NGROUPS];
  double integral_N[RT_NGROUPS];
  double integral_sigma_E[RT_NGROUPS][RT_NIONIZING_SPECIES];
  double integral_sigma_E_over_hnu[RT_NGROUPS][RT_NIONIZING_SPECIES];

  /* Prepare parameter struct for integration functions */
  /* -------------------------------------------------- */
  struct spectrum_integration_params integration_params = {/*species=*/0,
                                                           /*T=*/0.,
                                                           /*kB=*/0.,
                                                           /*h_planck=*/0.,
                                                           /*c=*/0.,
                                                           /*params=*/NULL};
  integration_params.T = T_blackbody;
  integration_params.kB = const_kboltz;
  integration_params.h_planck = const_planck_h;
  integration_params.c = const_speed_light_c;
  struct photoion_cs_parameters cs_params = init_photoion_cs_params_cgs();
  integration_params.cs_params = &cs_params;

  /* Compute integrals */
  /* ----------------- */

  double nu_start[RT_NGROUPS];
  double nu_stop[RT_NGROUPS];
  for (int group = 0; group < RT_NGROUPS; group++)
    nu_start[group] = frequency_bins[group];
  for (int group = 0; group < RT_NGROUPS - 1; group++)
    nu_stop[group] = frequency_bins[group + 1];

  /* don't start at exactly 0 to avoid unlucky divisions */
  if (nu_start[0] == 0.)
    nu_start[0] = nu_start[1] * TINY_NUMBER;
  if (RT_NGROUPS == 1) {
    /* If we only have one group, start integrating from the Hydrogen
     * ionization frequency, not from zero. */
    nu_start[0] = cs_params.E_ion[0] / const_planck_h;
  }
  /* Stop at 10 times the frequency of the blackbody peak */
  nu_stop[RT_NGROUPS - 1] =
      10. * blackbody_peak_frequency(T_blackbody, const_kboltz, const_planck_h);

  for (int group = 0; group < RT_NGROUPS; group++) {
    /* This is independent of species. */
    integral_E[group] = cross_sections_integrate_gsl(
        spectrum_integrand, nu_start[group], nu_stop[group],
        RT_INTEGRAL_NPOINTS, &integration_params);
    integral_N[group] = cross_sections_integrate_gsl(
        spectrum_over_hnu_integrand, nu_start[group], nu_stop[group],
        RT_INTEGRAL_NPOINTS, &integration_params);

    for (int species = 0; species < RT_NIONIZING_SPECIES; species++) {
      integration_params.species = species;
      integral_sigma_E[group][species] = cross_sections_integrate_gsl(
          spectrum_times_sigma_integrand, nu_start[group], nu_stop[group],
          RT_INTEGRAL_NPOINTS, &integration_params);
      integral_sigma_E_over_hnu[group][species] = cross_sections_integrate_gsl(
          spectrum_times_sigma_over_hnu_integrand, nu_start[group],
          nu_stop[group], RT_INTEGRAL_NPOINTS, &integration_params);
    }
  }

  /* Now compute the actual average cross sections */
  /* --------------------------------------------- */
  for (int group = 0; group < RT_NGROUPS; group++) {
    if (integral_N[group] > 0.) {
      mean_energy[group] = integral_E[group] / integral_N[group];
    } else {
      mean_energy[group] = 0.;
    }
    for (int spec = 0; spec < RT_NIONIZING_SPECIES; spec++) {
      if (integral_E[group] > 0.) {
        cse[group][spec] = integral_sigma_E[group][spec] / integral_E[group];
        csn[group][spec] =
            integral_sigma_E_over_hnu[group][spec] / integral_N[group];
      } else {
        /* No radiation = no interaction */
        cse[group][spec] = 0.;
        csn[group][spec] = 0.;
      }
    }
  }
}

#endif /* CROSS_SECTIONS_H */
