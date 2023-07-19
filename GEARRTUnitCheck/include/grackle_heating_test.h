#ifndef GRACKLE_HEATING_TEST_H
#define GRACKLE_HEATING_TEST_H

/* ---------------------------------------------
 * Reproduce the Iliev06 Test 0 part 3 example:
 * Heat the gas with radiation for 0.5 Myr, then
 * let it cool for 5 more Myrs
 * Except we vary the initial density here.
 * --------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <grackle.h>

#include "grackle_checks.h"
#include "ionization_equilibrium.h"
#include "mean_molecular_weight.h"
#include "photon_interaction_rates.h"

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
 * @param density gas density to use
 * @param fixed_radiation_density_field_cgs fixed energy density to use as
 *heating source, in erg/cm^3/s
 * @param name name of the test case. NO SPACES.
 * @param mass_units the internal mass units to use.
 * @param length_units the internal length units to use.
 * @param density_units the internal density units to use.
 * @param velocity_units the internal velocity units to use.
 * @param internal_energy_units the internal specific internal energy units
 * @param dump_results whether to write results to file
 * @param verbose print time step data to screen?
 **/
void run_grackle_heating_test(
    float density, double fixed_radiation_density_field_cgs[RT_NGROUPS],
    char *name, double mass_units, double length_units, double time_units,
    double density_units, double velocity_units, double internal_energy_units,
    int dump_results, int verbose) {

  /******************************************************************
   * Set up initial conditions and runtime parameters.
   *****************************************************************/

  message("Running grackle heating test for case %s", name);

  /* output file */
  char filename[200];
  sprintf(filename, "heating_test-%s.dat", name);
  FILE *fd = NULL;
  if (dump_results)
    fd = fopen(filename, "w");
  /* output frequency  in number of steps */
  int output_frequency_cool = 64; /* output frequency while cooling */
  int output_frequency_heat = 2;  /* output frequency while heating */

  /* Time integration variables */
  /* -------------------------- */
  /* max dt while cooling. in yr. Will be converted later */
  double dt_max_cool = 2000.;
  /* max dt while heating. in yr. Will be converted later */
  double dt_max_heat = 100.;
  double tinit = 1e-5; /* in yr; will be converted later */
  double tend = 5.5e6; /* in yr; will be converted later */

  double t = tinit * const_yr / time_units;          /* yr to code units */
  tend = tend * const_yr / time_units;               /* yr to code units */
  dt_max_cool = dt_max_cool * const_yr / time_units; /* yr to code units */
  dt_max_heat = dt_max_heat * const_yr / time_units; /* yr to code units */

  /* Set up initial conditions for gas cells */
  /* --------------------------------------- */
  double hydrogen_fraction_by_mass = 1.00;
  double gas_density = density;

  /* Derived quantities from ICs */
  /* --------------------------- */

  /* Assuming fully neutral hydrogen gas */
  double mu_init = 1.;
  double T = 100.; /* K */
  double internal_energy_cgs =
      T * const_kboltz / (const_mh * mu_init * (const_adiabatic_index - 1.));
  double internal_energy = internal_energy_cgs / internal_energy_units;

  /* define the hydrogen number density */
  /* use `gr_float` to use the same type of floats that grackle
   * is compiled in. Compile grackle with precision-32 if you want
   * floats, or precision-64 otherwise. */
  gr_float nH =
      hydrogen_fraction_by_mass * gas_density / (const_mh / mass_units);
  gr_float nHI, nHII, nHeI, nHeII, nHeIII, ne;

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
  gr_float species_densities[12] = {HI_density,
                                    HII_density,
                                    HeI_density,
                                    HeII_density,
                                    HeIII_density,
                                    e_density,
                                    TINY_NUMBER * gas_density,
                                    TINY_NUMBER * gas_density,
                                    TINY_NUMBER * gas_density,
                                    TINY_NUMBER * gas_density,
                                    TINY_NUMBER * gas_density,
                                    TINY_NUMBER * gas_density};

  /* Get photon cross sections and mean energies */
  /* ------------------------------------------- */

#if RT_NGROUPS == 3
  /* Just store the solutions, don't mind the computation for this test.
   * It's always done in cgs, so that shouldn't fail regardless of user
   * chosen units. */
  /* Cross sections are in cm^-2, mean photon energies in erg. */
  double cse_0[3] = {2.782672e-18, 0., 0.};
  double cse_1[3] = {5.043639e-19, 6.020192e-24, 0.};
  double cse_2[3] = {7.459092e-20, 6.515781e-25, 1.002143e-18};
  double *cse[RT_NGROUPS] = {cse_0, cse_1, cse_2};
  double csn_0[3] = {3.007890e-18, 0., 0.};
  double csn_1[3] = {5.688884e-19, 6.901561e-24, 0.};
  double csn_2[3] = {7.893315e-20, 6.929928e-25, 1.056205e-18};
  double *csn[RT_NGROUPS] = {csn_0, csn_1, csn_2};
  double mean_photon_energies[RT_NGROUPS] = {3.020e-11, 5.619e-11, 1.052e-10};
#elif RT_NGROUPS == 1
  double cse_0[3] = {1.096971e-18, 2.303147e-24, 1.298503e-19};
  double *cse[RT_NGROUPS] = {cse_0};
  double csn_0[3] = {1.630511e-18, 3.450498e-24, 6.171427e-20};
  double *csn[RT_NGROUPS] = {csn_0};
  double mean_photon_energies[RT_NGROUPS] = {4.744e-11};
#endif

  double radiation_energy_density_cgs[RT_NGROUPS];
  for (int g = 0; g < RT_NGROUPS; g++) {
    radiation_energy_density_cgs[g] = fixed_radiation_density_field_cgs[g];
  }

  gr_float ion_densities_cgs[6];
  /* I need them in cgs. */
  ion_densities_cgs[0] = species_densities[0] * density_units;
  ion_densities_cgs[1] = species_densities[1] * density_units;
  ion_densities_cgs[2] = species_densities[2] * density_units;
  ion_densities_cgs[3] = species_densities[3] * density_units;
  ion_densities_cgs[4] = species_densities[4] * density_units;
  ion_densities_cgs[5] = species_densities[5] * density_units;

  gr_float interaction_rates[5] = {0., 0., 0., 0., 0};

  get_interaction_rates(radiation_energy_density_cgs, ion_densities_cgs, cse,
                        csn, mean_photon_energies, time_units,
                        interaction_rates);

  /* Grackle behaviour setup */
  /* ----------------------- */

  int UVbackground = 0;         /* toogle on/off the UV background */
  int primordial_chemistry = 1; /* choose the chemical network */
  int use_radiative_cooling = 1;
  int use_radiative_transfer = 1;
  char *grackle_data_file = "";

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
  if (local_initialize_chemistry_parameters(&grackle_chemistry_data) == 0) {
    error("Error in set_default_chemistry_parameters.");
  }

  /* Set parameter values for chemistry. */
  setup_grackle_chemistry(&grackle_chemistry_data, primordial_chemistry,
                          UVbackground, grackle_data_file,
                          use_radiative_cooling, use_radiative_transfer,
                          hydrogen_fraction_by_mass);

  chemistry_data_storage grackle_chemistry_rates;

  if (local_initialize_chemistry_data(&grackle_chemistry_data,
                                      &grackle_chemistry_rates,
                                      &grackle_units_data) == 0) {
    error("Error in local_initialize_chemistry_data");
  }

  /* Gas Data */
  /* -------- */

  /* Create struct for storing grackle field data */
  grackle_field_data grackle_fields;
  setup_grackle_fields(&grackle_fields, species_densities, interaction_rates,
                       gas_density, internal_energy);

  /* Write headers */
  /* ------------- */

  if (verbose) {
    write_header(stdout);
    write_timestep(stdout, &grackle_fields, &grackle_units_data,
                   &grackle_chemistry_data, &grackle_chemistry_rates,
                   /*field_index=*/0, t, dt_max_heat, time_units, /*step=*/0);
  }

  /* write down what ICs you used into file */
  if (dump_results) {
    write_my_setup(fd, grackle_fields, &grackle_chemistry_data, mass_units,
                   length_units, velocity_units, dt_max_heat,
                   hydrogen_fraction_by_mass, gas_density, internal_energy);
    write_header(fd);
    write_timestep(fd, &grackle_fields, &grackle_units_data,
                   &grackle_chemistry_data, &grackle_chemistry_rates,
                   /*field_index=*/0, t, dt_max_heat, time_units, /*step=*/0);
  }

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/

  int step = 0;
  int completion = 0;
  int completion_fractions = 5;
  while (t < tend) {

    /* Set up radiation fields, and compute the resulting interaction
     * rates depending on the simulation time. */
    gr_float iact_rates[5] = {0., 0., 0., 0., 0.};
    double dt_max = dt_max_cool;
    int output_frequency_to_use = output_frequency_cool;
    if (t / const_yr * time_units < 0.5e6) {
      /* below 0.5 Myr, we heat. */
      iact_rates[0] = interaction_rates[0];
      iact_rates[1] = interaction_rates[1];
      iact_rates[2] = interaction_rates[2];
      iact_rates[3] = interaction_rates[3];
      iact_rates[4] = interaction_rates[4];
      dt_max = dt_max_heat;
      output_frequency_to_use = output_frequency_heat;
    }

    /* Tell grackle about the rates */
    for (int i = 0; i < FIELD_SIZE; i++) {
      grackle_fields.RT_heating_rate[i] = iact_rates[0];
      grackle_fields.RT_HI_ionization_rate[i] = iact_rates[1];
      grackle_fields.RT_HeI_ionization_rate[i] = iact_rates[2];
      grackle_fields.RT_HeII_ionization_rate[i] = iact_rates[3];
      grackle_fields.RT_H2_dissociation_rate[i] = iact_rates[4];
    }

    /* Get cooling time */
    gr_float tchem_time;
    if (local_calculate_cooling_time(
            &grackle_chemistry_data, &grackle_chemistry_rates,
            &grackle_units_data, &grackle_fields, &tchem_time) == 0)
      error("Error in calculate_cooling_time.");
    double dt = fmin(fabs(tchem_time), dt_max);

    t += dt;
    step += 1;

    if (local_solve_chemistry(&grackle_chemistry_data, &grackle_chemistry_rates,
                              &grackle_units_data, &grackle_fields, dt) == 0) {
      error("Error in solve_chemistry.");
    }

    grackle_checks_density_sum(density, &grackle_fields);
    grackle_checks_ion_sum(&grackle_fields, mass_units);

    if (verbose) {
      write_timestep(stdout, &grackle_fields, &grackle_units_data,
                     &grackle_chemistry_data, &grackle_chemistry_rates,
                     /*field_index=*/0, t, dt, time_units, step);
    } else {
      if (t / tend >
          ((double)(completion + 1) / (double)completion_fractions)) {
        message("Completed %.1f %%",
                (float)(completion + 1) / (float)completion_fractions * 100.f);
        completion++;
      }
    }

    if (step % output_frequency_to_use == 0 && dump_results)
      write_timestep(fd, &grackle_fields, &grackle_units_data,
                     &grackle_chemistry_data, &grackle_chemistry_rates,
                     /*field_index=*/0, t, dt, time_units, step);
  }

  /* Cleanup */
  if (dump_results)
    fclose(fd);
  clean_up_fields(&grackle_fields);
  local_free_chemistry_data(&grackle_chemistry_data, &grackle_chemistry_rates);

  return;
}
#endif /* GRACKLE_HEATING_TEST_H */
