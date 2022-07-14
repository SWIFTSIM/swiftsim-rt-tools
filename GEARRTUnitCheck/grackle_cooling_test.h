#ifndef GRACKLE_COOLING_TEST_H
#define GRACKLE_COOLING_TEST_H

/* ----------------------------------------------------------
 * In this example, we start with high internal energies and 
 * a fully ionized gas, and just let it cool without any RT.
 * --------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <grackle.h>

#include "ionization_equilibrium.h"
#include "mean_molecular_weight.h"
#include "error.h"
#include "grackle_checks.h"
/* NOTE: don't include my_grackle_utils here to make sure you've
 * got the correct FIELD_SIZE etc definitions. */
#ifndef MY_GRACKLE_UTILS_H
#error This file needs my_grackle_utils.h to be included already, \
  but not from within the file itself.
#endif

/**
 * @brief Run a test with an initially high temperature gas
 * that cools down.
 * For the units, I use the following convention:
 * quantity * units = quantity_in_cgs
 *
 * @param denstiy gas density to use
 * @param name name of the test case. NO SPACES.
 * @param mass_units the internal mass units to use. 
 * @param length_units the internal length units to use. 
 * @param density_units the internal density units to use. 
 * @param velocity_units the internal velocity units to use. 
 * @param internal_energy_units the internal specific internal energy units to use. 
 * @param verbose print time step data to screen?
 **/
void run_grackle_cooling_test(
    double density,
    char *name,
    double mass_units,
    double length_units,
    double time_units,
    double density_units,
    double velocity_units,
    double internal_energy_units,
    int verbose
    ) {

  /******************************************************************
   * Set up initial conditions and runtime parameters.
   *****************************************************************/

  message("Running grackle cooling test for case %s", name);

  /* output file */
  char filename[200];
  sprintf(filename, "cooling_test-%s.dat", name);
  message("filename %s", filename);
  FILE *fd = fopen(filename, "w");
  /* output frequency  in number of steps */
  int output_frequency = 8;

  /* Time integration variables */
  /* -------------------------- */
  double dt = 1e3;    /* in yr. will be converted later */
  double tinit = 1.;  /* in yr; will be converted later */
  double tend = 1e7;  /* in yr; will be converted later */

  double t = tinit * const_yr / time_units;   /* yr to code units */
  dt = dt * const_yr / time_units;   /* yr to code units */
  tend = tend * const_yr / time_units; /* yr to code units */

  /* Set up initial conditions for gas cells */
  /* --------------------------------------- */
  double hydrogen_fraction_by_mass = 0.76;
  /* Use solution of swift's output. This is in internal units already. */
  double gas_density = density;

  /* Derived quantities from ICs */
  /* --------------------------- */

  /* Assuming fully ionized gas */
  double mu_init = mean_molecular_weight_from_mass_fractions(
      0., hydrogen_fraction_by_mass, 0., 0., (1. - hydrogen_fraction_by_mass));

  double T = 1e4; /* K */
  double internal_energy_cgs = const_kboltz * T / ((const_adiabatic_index - 1.) * mu_init * const_mh);
  double internal_energy = internal_energy_cgs / internal_energy_units;

  /* define the hydrogen number density */
  /* use `gr_float` to use the same type of floats that grackle
   * is compiled in. Compile grackle with precision-32 if you want
   * floats, or precision-64 otherwise. */
  gr_float nH =
      hydrogen_fraction_by_mass * gas_density / (const_mh / mass_units);
  gr_float nHI;
  gr_float nHII;
  gr_float nHeI;
  gr_float nHeII;
  gr_float nHeIII;
  gr_float ne;

  /* get densities of primordial spicies assuming ionization equilibrium */
  ionization_equilibrium_calculate_densities(T, nH, hydrogen_fraction_by_mass,
                                             &nHI, &nHII, &nHeI, &nHeII,
                                             &nHeIII, &ne);

  gr_float HI_density = nHI * (const_mh / mass_units);
  gr_float HII_density = nHII * (const_mh / mass_units);
  gr_float HeI_density = nHeI * (4 * const_mh / mass_units);
  gr_float HeII_density = nHeII * (4 * const_mh / mass_units);
  gr_float HeIII_density = nHeIII * (4 * const_mh / mass_units);
  /* !! this is the convention adopted by Grackle
   * !! e_density is the electron number density multiplied by proton mass,
   * !! or electron mass density * nH / ne */
  gr_float e_density = ne * (const_mh / mass_units);

  /* Store them all in a single array for simplicity. */
  gr_float species_densities[12] = {
      HI_density, HII_density, HeI_density, HeII_density, HeIII_density,
      e_density,  TINY_NUMBER, TINY_NUMBER, TINY_NUMBER, TINY_NUMBER,
      TINY_NUMBER, TINY_NUMBER};

  /* Grackle behaviour setup */
  /* ----------------------- */

  int UVbackground = 0;         /* toogle on/off the UV background */
  int primordial_chemistry = 1; /* choose the chemical network */
  int use_radiative_cooling = 1;
  int use_radiative_transfer = 0;
  char *grackle_data_file = "";

  gr_float RT_HI_ionization_rate = 0.;
  gr_float RT_HeI_ionization_rate = 0.;
  gr_float RT_HeII_ionization_rate = 0.;
  gr_float RT_H2_dissociation_rate = 0.;
  gr_float RT_heating_rate = 0.;

  /* Store them all in a single array for simplicity. */
  gr_float interaction_rates[5] = {
      RT_HI_ionization_rate, RT_HeI_ionization_rate, RT_HeII_ionization_rate,
      RT_H2_dissociation_rate, RT_heating_rate};

  /*********************************************************************
   * Set up gracke data and fields.
   **********************************************************************/

  /* Units  */
  /* ------ */

  /* First, set up the units system. We assume cgs
   * These are conversions from code units to cgs. */
  code_units grackle_units_data;
  setup_grackle_units(&grackle_units_data, density_units, length_units,
                      time_units);

  /* Chemistry Parameters */
  /* -------------------- */

  /* Second, create a chemistry object for parameters.  This needs to be a
   * pointer. */
  chemistry_data grackle_chemistry_data;
  if (set_default_chemistry_parameters(&grackle_chemistry_data) == 0) {
    fprintf(stderr, "Errir in set_default_chemistry_parameters");
    abort();
  }

  /* Set parameter values for chemistry. */
  setup_grackle_chemistry(&grackle_chemistry_data, primordial_chemistry,
                          UVbackground, grackle_data_file,
                          use_radiative_cooling, use_radiative_transfer,
                          hydrogen_fraction_by_mass);

  /* Initialize the chemistry_data_storage object to be able to use local
   * functions */
  if (initialize_chemistry_data(&grackle_units_data) == 0) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    abort();
  }

  /* Gas Data */
  /* -------- */

  /* Create struct for storing grackle field data */
  grackle_field_data grackle_fields;
  setup_grackle_fields(&grackle_fields, species_densities, interaction_rates,
                       gas_density, internal_energy);

  /* Write headers */
  /* ------------- */

  /* First to stdout */

  /* Now into a file as well. */
  /* also write down what ICs you used into file */
  write_my_setup(fd, grackle_fields, grackle_chemistry_data, mass_units,
                 length_units, velocity_units, dt, hydrogen_fraction_by_mass,
                 gas_density, internal_energy);
  write_header(fd);
  write_timestep(fd, &grackle_fields, &grackle_units_data,
                 &grackle_chemistry_data, /*field_index=*/0, t, dt, time_units,
                 /*step=*/0);

  if (verbose){
    /* Print to screen as well */
    write_header(stdout);
    write_timestep(stdout, &grackle_fields, &grackle_units_data, &grackle_chemistry_data, /*field_index=*/0, t, dt, time_units, /*step=*/0);
  }

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/

  int step = 0;
  while (t < tend) {

    t += dt;
    step += 1;

    /* if (solve_chemistry(&grackle_units_data, &grackle_fields, dt) == 0) { */
    if (local_solve_chemistry(&grackle_chemistry_data, &grackle_rates,
                              &grackle_units_data, &grackle_fields, dt) == 0) {
      fprintf(stderr, "Error in solve_chemistry.\n");
      abort();
    }

    grackle_checks_density_sum(density, &grackle_fields);

    if (verbose)
      write_timestep(stdout, &grackle_fields, &grackle_units_data,
                     &grackle_chemistry_data, /*field_index=*/0, t, dt,
                     time_units, step);
    if (step % output_frequency == 0)
      write_timestep(fd, &grackle_fields, &grackle_units_data,
                     &grackle_chemistry_data, /*field_index=*/0, t, dt,
                     time_units, step);
  }

  /* Clean up after yourself */

  fclose(fd);
  clean_up_fields(&grackle_fields);

  return;
}
#endif
