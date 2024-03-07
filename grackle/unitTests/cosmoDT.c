/* ---------------------------------------------
 * Compute total time between a_begin and a_end
 * as a series of intervals computing dt.
 * Then make sure it adds up correctly.
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

  /* Cosmology                  */
  /* -------------------------- */
  double a_begin = 1e-10;
  double a_end = 1.;
  int nsteps = 1000000;

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

  double a_table[COSMO_TABLE_ELEMENTS];
  double t_table[COSMO_TABLE_ELEMENTS];

  cosmo_get_tables(a_table, t_table, &cosmology, a_begin, a_end);

  const double t_tot = cosmo_get_dt(a_begin, a_end, a_begin, a_end, t_table);

  double a;
  double a_next;

  const double da = (a_end - a_begin) / nsteps;
  const double log_a_begin = log(a_begin);
  const double log_a_end = log(a_end);
  const double dlog_a = (log_a_end - log_a_begin) / nsteps;

  double tolerance = 1.e-4;
  if (10. / (double)nsteps > tolerance)
    tolerance = 10. / (double)nsteps;

  /* Use linear steps */
  /* ---------------- */
  double t_linear = 0.;

  a = a_begin;
  int steps = 0;
  while (a < a_end) {
    steps++;
    a_next = a + da;
    double dt = cosmo_get_dt(a, a_next, a_begin, a_end, t_table);

    t_linear += dt;
    a = a_next;
  }

  printf("LINEAR:      Expect=%12.6g Got=%12.6g Ratio=%12.6g Steps=%6i\n",
         t_tot, t_linear, t_linear / t_tot, steps);

  double dev = fabs(1. - t_linear / t_tot);
  if (dev > tolerance)
    error("Deviation too large: %g", dev);

  /* Use log steps */
  /* ---------------- */
  double t_log = 0.;

  double log_a = log_a_begin;
  a = a_begin;
  steps = 0;
  while (a < a_end) {
    steps++;

    double log_a_next = log_a + dlog_a;
    a_next = exp(log_a_next);
    double dt = cosmo_get_dt(a, a_next, a_begin, a_end, t_table);

    t_log += dt;
    a = a_next;
    log_a = log_a_next;
  }

  printf("LOGARITHMIC: Expect=%12.6g Got=%12.6g Ratio=%12.6g Steps=%6i\n",
         t_tot, t_log, t_log / t_tot, steps);

  dev = fabs(1. - t_log / t_tot);
  if (dev > tolerance)
    error("Deviation too large: %g", dev);

  printf("Done.\n");

  return 0;
}
