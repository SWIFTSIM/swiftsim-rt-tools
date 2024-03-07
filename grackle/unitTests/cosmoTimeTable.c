/* ---------------------------------------------
 * Compute the a(t) and t(a) time tables and
 * write them out.
 * --------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "cosmology.h"

int main() {

  /* Define units : use the same as internal units for swift */
  /* ------------------------------------------------------- */
  /* double mass_units = 1.99848e43; */
  double length_units = 3.08567758e21;
  double velocity_units = 1e5;

  /* double density_units = mass_units / length_units / length_units /
   * length_units; */
  double time_units = length_units / velocity_units;

  /* time_units = 2.; */

  /* Cosmology                  */
  /* -------------------------- */
  double a_begin = 1e-10;
  double a_end = 1.;

  struct cosmology cosmology;
  /* Planck13 (EAGLE flavour) */

  cosmology.Omega_cdm = 0.2587; /* Dark matter density parameter*/
  cosmology.Omega_b = 0.04825;  /* baryon density parameter*/
  cosmology.Omega_l = 0.693;    /* Dark Energy density parameter */
  cosmology.Omega_k = 0.;       /* Radiation density parameter */
  cosmology.Omega_r = 0.;       /* Radiation density parameter */
  cosmology.Omega_nu = 0.;      /* Neutrino density parameter */
  cosmology.w_0 = -1.0; /* Dark-energy equation-of-state parameter at z=0. */
  cosmology.w_a =
      0.; /* Dark-energy equation-of-state time evolution parameter. */
  cosmology.H_0 = 67.77; /* Hubble constant at z=0 in km/s/Mpc */

  cosmo_convert_H0_to_internal_units(&cosmology, time_units);

  /* Compute a(t) and t(a) tables for interpolation */
  /* ---------------------------------------------- */

  double *a_table = malloc(sizeof(double) * COSMO_TABLE_ELEMENTS);
  double *t_table = malloc(sizeof(double) * COSMO_TABLE_ELEMENTS);

  cosmo_get_tables(a_table, t_table, &cosmology, a_begin, a_end);

  /* output file */
  FILE *fd = fopen("cosmoTimeTable.dat", "w");
  fprintf(fd, "#%11s, %12s\n", "t [yr]", "a [1]");
  for (int i = 0; i < COSMO_TABLE_ELEMENTS; i++) {
    double t = t_table[i] / const_yr * time_units;
    double a = a_table[i];
    fprintf(fd, "%12.6g, %12.6g\n", t, a);
  }
  fclose(fd);

  free(a_table);
  free(t_table);

  printf("Done.\n");

  return 0;
}
