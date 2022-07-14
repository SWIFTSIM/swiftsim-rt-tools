/* ---------------------------------------------
 * In this example, we start with high internal
 * energies and a fully ionized gas, and just
 * let it cool without any RT.
 * --------------------------------------------- */

/* define these before including local headers like my_grackle_utils.h */
#define FIELD_SIZE 1
#define GRIDDIM 1

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <grackle.h>

#include "constants.h"
#include "cross_sections.h"
#include "ionization_equilibrium.h"
#include "mean_molecular_weight.h"
#include "my_grackle_utils.h"
#include "photon_interaction_rates.h"

int main() {

  /******************************************************************
   * Set up initial conditions and runtime parameters.
   *****************************************************************/

  /* Print some extra data to screen? */
  int verbose = 1;
  grackle_verbose = 0;
  /* output file */
  FILE *fd = fopen("out.dat", "w");
  /* output frequency  in number of steps */
  int output_frequency_cool = 64; /* output frequency while cooling */
  int output_frequency_heat = 2;  /* output frequency while heating */

  /* Define units  */
  /* --------------*/
  /* Stick with cgs for now. */
  /* double mass_units = const_mh; */
  double mass_units = 1.67262171e-24;
  /* double length_units = 1e3; */
  /* double velocity_units = 1e-6; */
  /* double mass_units = 1.e-10; */
  double length_units = 1.;
  double velocity_units = 1.;

  double density_units =
      mass_units / length_units / length_units / length_units;
  double time_units = length_units / velocity_units;

  /* Time integration variables */
  /* -------------------------- */
  double dt_max_cool =
      1000.; /* max dt while cooling. in yr. Will be converted later */
  double dt_max_heat =
      2.; /* max dt while heating. in yr. Will be converted later */
  double tinit = 1e-5; /* in yr; will be converted later */
  double tend = 5.5e6; /* in yr; will be converted later */
  /* Convert times to internal units. */
  double t = tinit * const_yr / time_units;          /* yr to code units */
  tend = tend * const_yr / time_units;               /* yr to code units */
  dt_max_cool = dt_max_cool * const_yr / time_units; /* yr to code units */
  dt_max_heat = dt_max_heat * const_yr / time_units; /* yr to code units */

  /* Set up initial conditions for gas cells */
  /* --------------------------------------- */
  double hydrogen_fraction_by_mass = 1.00;
  double gas_density = const_mh; /* in cgs, will be converted later.
                                    Corresponds to number density 1 cm^-3 */
  double T = 100.;               /* K */

  gas_density /= density_units;

  /* Initial conditions for radiation */
  /* See README for details */
  double T_blackbody = 1e5; /* K */
#if RT_NGROUPS == 4
  double frequency_bins_Hz[4] = {0., 3.288e15, 5.945e15, 13.157e15}; /* Hz */
  double fixed_luminosity_cgs[4] = {0., 1.350e+01, 2.779e+01,
                                    6.152e+00}; /* erg / cm^2 / s */
#elif RT_NGROUPS == 1
  double frequency_bins_Hz[1] = {3.288e15};     /* Hz */
  double fixed_luminosity_cgs[1] = {4.774e+01}; /* erg / cm^2 / s */
#else
  printf("You need to set up the correct frequency bins and luminosities "
         "for " RT_NGROUPS " groups used");
  return EXIT_FAILURE;
#endif
  double fixed_radiation_density_field_cgs[4] = {0., 0., 0., 0.};
  for (int g = 0; g < RT_NGROUPS; g++) {
    fixed_radiation_density_field_cgs[g] =
        fixed_luminosity_cgs[g] / const_speed_light_c;
  }

  /* Derived quantities from ICs */
  /* --------------------------- */

  /* Assuming fully neutral gas */
  double mu_init = 1.;
  double internal_energy =
      T * const_kboltz / (const_mh * mu_init * (const_adiabatic_index - 1.));

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
   * !! e_density is the electron number density multiplied by proton mass */
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
  /* Note that the result is always in cgs. */
  double **cse = malloc(RT_NGROUPS * sizeof(double *));
  double **csn = malloc(RT_NGROUPS * sizeof(double *));
  double mean_photon_energies[RT_NGROUPS];
  for (int group = 0; group < RT_NGROUPS; group++) {
    cse[group] = malloc(RT_NIONIZING_SPECIES * sizeof(double));
    csn[group] = malloc(RT_NIONIZING_SPECIES * sizeof(double));
    mean_photon_energies[group] = 0.;
  }

  get_cross_sections(T_blackbody, frequency_bins_Hz, cse, csn,
                     mean_photon_energies);

  const double eV_per_erg = 6.242e+11;
  for (int g = 0; g < RT_NGROUPS; g++) {
    printf("group %d cse=%.3e csn=%.3e mean=%.3e erg %.3e eV\n", g, cse[g][0],
           csn[g][0], mean_photon_energies[g],
           mean_photon_energies[g] * eV_per_erg);
  }

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
  printf("Mladen's rates %12.3e %12.3e %12.3e %12.3e %12.3e\n",
         interaction_rates[0], interaction_rates[1], interaction_rates[2],
         interaction_rates[3], interaction_rates[4]);

  /* Grackle behaviour setup */
  /* ----------------------- */

  int UVbackground = 0;         /* toogle on/off the UV background */
  int primordial_chemistry = 1; /* choose the chemical network */
  int use_radiative_cooling = 1;
  int use_radiative_transfer = 1; // YR: was 0
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
  if (set_default_chemistry_parameters(&grackle_chemistry_data) == 0) {
    fprintf(stderr, "Error in set_default_chemistry_parameters");
    return EXIT_FAILURE;
  }

  /* Set parameter values for chemistry. */
  setup_grackle_chemistry(&grackle_chemistry_data, primordial_chemistry,
                          UVbackground, grackle_data_file,
                          use_radiative_cooling, use_radiative_transfer,
                          hydrogen_fraction_by_mass);

  if (initialize_chemistry_data(&grackle_units_data) == 0) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return EXIT_FAILURE;
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
    printf("%15s%15s%15s%15s%15s%15s%15s%15s\n",
           "Initial setup: ", "Temperature", "nHI", "nHII", "nHeI", "nHeII",
           "nHeIII", "ne");
    printf("%15s%15g%15g%15g%15g%15g%15g%15g\n\n", "Initial setup: ", T, nHI,
           nHII, nHeI, nHeII, nHeIII, ne);
  }

  write_header(stdout);
  write_timestep(stdout, &grackle_fields, &grackle_units_data,
                 &grackle_chemistry_data, /*field_index=*/0, t, dt_max_heat,
                 time_units, /*step=*/0);

  /* write down what ICs you used into file */
  /* TODO: dt is wrong here */
  write_my_setup(fd, grackle_fields, grackle_chemistry_data, mass_units,
                 length_units, velocity_units, dt_max_heat,
                 hydrogen_fraction_by_mass, gas_density, internal_energy);
  write_header(fd);
  write_timestep(fd, &grackle_fields, &grackle_units_data,
                 &grackle_chemistry_data, /*field_index=*/0, t, dt_max_heat,
                 time_units, /*step=*/0);

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/

  int step = 0;
  while (t < tend) {

    /* Set up radiation fields, and compute the resulting interaction
     * rates depending on the simulation time. */

    gr_float iact_rates[5] = {0., 0., 0., 0., 0.};
    double dt = dt_max_cool;
    int output_frequency_to_use = output_frequency_cool;
    if (t / const_yr * time_units < 0.5e6) {
      /* below 0.5 Myr, we heat. */
      iact_rates[0] = interaction_rates[0];
      iact_rates[1] = interaction_rates[1];
      iact_rates[2] = interaction_rates[2];
      iact_rates[3] = interaction_rates[3];
      iact_rates[4] = interaction_rates[4];
      dt = dt_max_heat;
      output_frequency_to_use = output_frequency_heat;
    }

    /* printf("heating\n"); */

    //      double radiation_energy_density[RT_NGROUPS];
    //      radiation_energy_density[0] = fixed_radiation_density_field[0];
    //      radiation_energy_density[1] = fixed_radiation_density_field[1];
    //      radiation_energy_density[2] = fixed_radiation_density_field[2];
    //      radiation_energy_density[3] = fixed_radiation_density_field[3];

    //      gr_float species_densities[6];
    //      species_densities[0] = grackle_fields.HI_density[0];
    //      species_densities[1] = grackle_fields.HII_density[0];
    //      species_densities[2] = grackle_fields.HeI_density[0];
    //      species_densities[3] = grackle_fields.HeII_density[0];
    //      species_densities[4] = grackle_fields.HeIII_density[0];
    //      species_densities[5] = grackle_fields.e_density[0];

    //     get_interaction_rates(radiation_energy_density,
    //                           species_densities,
    //                           cse, csn, mean_photon_energies,
    //                           interaction_rates);
    /* } */

    for (int i = 0; i < FIELD_SIZE; i++) {
      /* grackle_fields.RT_heating_rate[i] = interaction_rates_yves[0]; */
      /* grackle_fields.RT_HI_ionization_rate[i] = interaction_rates_yves[1]; */
      /* grackle_fields.RT_HeI_ionization_rate[i] = interaction_rates_yves[2];
       */
      /* grackle_fields.RT_HeII_ionization_rate[i] = interaction_rates_yves[3];
       */
      /* grackle_fields.RT_H2_dissociation_rate[i] = interaction_rates_yves[4];
       */
      grackle_fields.RT_heating_rate[i] = iact_rates[0];
      grackle_fields.RT_HI_ionization_rate[i] = iact_rates[1];
      grackle_fields.RT_HeI_ionization_rate[i] = iact_rates[2];
      grackle_fields.RT_HeII_ionization_rate[i] = iact_rates[3];
      grackle_fields.RT_H2_dissociation_rate[i] = iact_rates[4];
    }

    t += dt;
    step += 1;

    /* if (solve_chemistry(&grackle_units_data, &grackle_fields, dt) == 0) { */
    if (local_solve_chemistry(&grackle_chemistry_data, &grackle_rates,
                              &grackle_units_data, &grackle_fields, dt) == 0) {
      /* if (local_solve_chemistry(&grackle_chemistry_data,
       * grackle_chemistry_rates, &grackle_units_data, &grackle_fields, dt) ==
       * 0) { */
      fprintf(stderr, "Error in solve_chemistry.\n");
      return EXIT_FAILURE;
    }

    write_timestep(stdout, &grackle_fields, &grackle_units_data,
                   &grackle_chemistry_data, /*field_index=*/0, t, dt,
                   time_units, step);

    if (step % output_frequency_to_use == 0)
      write_timestep(fd, &grackle_fields, &grackle_units_data,
                     &grackle_chemistry_data, /*field_index=*/0, t, dt,
                     time_units, step);
  }

  /* Cleanup */
  fclose(fd);
  clean_up_fields(&grackle_fields);
  _free_chemistry_data(&grackle_chemistry_data, &grackle_rates);
  for (int g = 0; g < RT_NGROUPS; g++) {
    free(cse[g]);
    free(csn[g]);
  }
  free(cse);
  free(csn);

  return EXIT_SUCCESS;
}