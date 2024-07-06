#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "cosmology.h"
#include "redshift.h"

int main() {
  /* Temperature */
  const double BB_temp = 1e5;

  /* Arrays to keep track of different integrals and frequencies */
  double photon_energy[RT_NGROUPS];
  double integral_E[RT_NGROUPS];
  double frequency_bins[RT_NGROUPS];

  /* Frequency @ BB peak */
  double BB_peak =
      blackbody_peak_frequency(BB_temp, const_kboltz, const_planck_h);
  double delta_nu = 10 * BB_peak / RT_NGROUPS;

  for (int g = 0; g < RT_NGROUPS; g++) {
    frequency_bins[g] = delta_nu * g;
  }

  double nu_start[RT_NGROUPS];
  double nu_stop[RT_NGROUPS];
  for (int group = 0; group < RT_NGROUPS; group++) {
    nu_start[group] = frequency_bins[group];
  }
  for (int group = 0; group < RT_NGROUPS - 1; group++) {
    nu_stop[group] = frequency_bins[group + 1];
  }
  double nu_stop_final =
      10 * blackbody_peak_frequency(BB_temp, const_kboltz, const_planck_h);
  nu_stop[RT_NGROUPS - 1] = nu_stop_final;

  struct integration_params params = {BB_temp};

  double total_energy = 0;
  for (int g = 0; g < RT_NGROUPS; g++) {
    integral_E[g] =
        redshift_integrate_gsl(rt_get_energy_density, nu_start[g], nu_stop[g],
                               N_INTEGRAL_POINTS, &params);

    photon_energy[g] = 1. * integral_E[g];

    printf("[group %d] : frequency start %.10e\n", g + 1, nu_start[g]);
    printf("[group %d] : frequency stop %.10e\n", g + 1, nu_stop[g]);
    printf("[group %d] : photon_energy %.10e\n", g + 1, photon_energy[g]);
    printf("[group %d] : integral_E %.10e\n\n", g + 1, integral_E[g]);
    total_energy += photon_energy[g];
  }

  printf("Total energy: %.10e\n", total_energy);

  FILE *energies_file = fopen("photon_energies", "wb");
  FILE *bins_file = fopen("frequency_bins", "wb");

  fwrite(photon_energy, sizeof(double), RT_NGROUPS, energies_file);
  fwrite(frequency_bins, sizeof(double), RT_NGROUPS, bins_file);
  fclose(energies_file);
  fclose(bins_file);
  return EXIT_SUCCESS;
}
