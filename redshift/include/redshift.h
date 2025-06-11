#include "blackbody.h"
#include "constants.h"
#include <gsl/gsl_integration.h>

#ifndef N_INTEGRAL_POINTS
#define N_INTEGRAL_POINTS 100000
#endif
struct integration_params {
  double BB_temp;
};

void write_cosmo_header(FILE *fd) {
  fprintf(fd, "#%8s %12s %15s %15s %15s\n", "Step", "a", "PhotonEnergy [IU]",
          "Volume[IU]", "BB temp [K]");
}

void write_cosmo_timestep(FILE *fd, int step, double a, double PhotonEnergy,
                          double volume, double BB_temperature,
                          int shifting_BB) {
  fprintf(fd, "%9d %12.6f %15.3e %15.3e %15.3e %d\n", step, a, PhotonEnergy,
          volume, BB_temperature, shifting_BB);
}

double get_spectrum_derivative(const double nu, const double T) {
  return blackbody_spectrum_intensity_first_derivative(
      nu, T, const_kboltz, const_planck_h, const_speed_light_c);
}

double rt_get_spectrum(double nu, void *params) {
  struct integration_params *p = (struct integration_params *)params;
  double T = p->BB_temp;
  return blackbody_spectrum_intensity(nu, T, const_kboltz, const_planck_h,
                                      const_speed_light_c);
}

double rt_get_energy_density(double nu, void *params) {
  struct integration_params *p = (struct integration_params *)params;
  double T = p->BB_temp;
  return blackbody_energy_density(nu, T, const_kboltz, const_planck_h,
                                  const_speed_light_c);
}

double spectrum_over_hnu_integrand(double nu, void *params) {
  struct integration_params *p = (struct integration_params *)params;
  const double E = nu * const_planck_h;
  const double J = rt_get_spectrum(nu, p);
  if (E > 0.) {
    return J / E;
  } else {
    return 0.;
  }
}

double spectrum_derivative_times_nu_integrand(double nu, void *params) {
  struct integration_params *p = (struct integration_params *)params;
  double T = p->BB_temp;
  const double J = get_spectrum_derivative(nu, T);

  return nu * J;
}

double redshift_integrate_gsl(double (*function)(double, void *),
                              double nu_start, double nu_stop, int npoints,
                              struct integration_params *params) {
  gsl_function F;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(npoints);
  double result, error;

  F.function = function;
  F.params = (void *)params;

  gsl_integration_qags(&F, nu_start, nu_stop, /*epsabs=*/0, /*epsrel=*/1e-7,
                       npoints, w, &result, &error);

  gsl_integration_workspace_free(w);

  return result;
}

void calculate_spectrum_energy_density(double *av_redshift_energy,
                                       double *nu_bins, double T_bb) {
  double nu_start[RT_NGROUPS];
  double nu_stop[RT_NGROUPS];

  struct integration_params params = {T_bb};

  for (int group = 0; group < RT_NGROUPS; group++) {
    nu_start[group] = nu_bins[group];
  }
  for (int group = 0; group < RT_NGROUPS - 1; group++) {
    nu_stop[group] = nu_bins[group + 1];
  }

  double nu_stop_final =
      10 * blackbody_peak_frequency(T_bb, const_kboltz, const_planck_h);
  nu_stop[RT_NGROUPS - 1] = nu_stop_final;

  for (int g = 0; g < RT_NGROUPS; g++) {
    av_redshift_energy[g] =
        redshift_integrate_gsl(rt_get_energy_density, nu_start[g], nu_stop[g],
                               N_INTEGRAL_POINTS, &params);
  }
}

void calculate_redshift_average_energy(double *av_redshift_energy,
                                       double *nu_bins, double T_bb) {
  double nu_start[RT_NGROUPS];
  double nu_stop[RT_NGROUPS];
  double integral_E[RT_NGROUPS];
  double integral_nudIdnu[RT_NGROUPS];

  struct integration_params params = {T_bb};

  /* Set energy bins for integration */
  for (int group = 0; group < RT_NGROUPS; group++) {
    nu_start[group] = nu_bins[group];
  }
  for (int group = 0; group < RT_NGROUPS - 1; group++) {
    nu_stop[group] = nu_bins[group + 1];
  }
  double nu_stop_final =
      10 * blackbody_peak_frequency(T_bb, const_kboltz, const_planck_h);
  nu_stop[RT_NGROUPS - 1] = nu_stop_final;

  for (int g = 0; g < RT_NGROUPS; g++) {
    integral_E[g] = redshift_integrate_gsl(
        rt_get_spectrum, nu_start[g], nu_stop[g], N_INTEGRAL_POINTS, &params);
    integral_nudIdnu[g] = redshift_integrate_gsl(
        spectrum_derivative_times_nu_integrand, nu_start[g], nu_stop[g],
        N_INTEGRAL_POINTS, &params);

    if (integral_E[g] != 0.) {
      av_redshift_energy[g] = integral_nudIdnu[g] / integral_E[g];
    } else {
      av_redshift_energy[g] = 0.;
    }
  }
}

/* void calculate_boundary_terms(double *boundary_terms, double* nu_bins, double
  T_bb) { double integral_E[RT_NGROUPS];

  double nu_start[RT_NGROUPS];
  double nu_stop[RT_NGROUPS];
  struct integration_params params = {T_bb};

  for (int group = 0; group < RT_NGROUPS; group++) {
    nu_start[group] = nu_bins[group];
  }
  for (int group = 0; group < RT_NGROUPS - 1; group++) {
    nu_stop[group] = nu_bins[group + 1];
  }

  double nu_stop_final = 10 * blackbody_peak_frequency(T_bb, const_kboltz,
  const_planck_h); nu_stop[RT_NGROUPS-1] = nu_stop_final;

  for (int g=0; g < RT_NGROUPS; g++) {
    integral_nuE[g] =
  redshift_integrate_gsl(spectrum_derivative_times_nu_integrand, nu_start[g],
  nu_stop[g], N_INTEGRAL_POINTS, &params);
  } */
