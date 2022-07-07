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
  int output_frequency = 1;

  /* Define units  */
  /* --------------*/
  /* Stick with cgs for now. */
  double mass_units = 1.;
  double length_units = 1.;
  double velocity_units = 1.;

  double density_units =
      mass_units / length_units / length_units / length_units;
  double time_units = length_units / velocity_units;

  /* Time integration variables */
  /* -------------------------- */
  double t = 0.;
  double dt_max = 5293212890.625; /* ~168yrs in internal units, 32768 time steps. */
  /* double dt_max = 330825805.6640625; [> 32768 * 16 <] */
  /* double dt = 661651611.328125; [> 32768 * 8 <] */
  double tinit = 1e-5;       /* in yr; will be converted later */
  double tend = 5.5e6;      /* in yr; will be converted later */
  t = tinit * const_yr / time_units;   /* yr to code units */
  tend = tend * const_yr / time_units; /* yr to code units */

  /* Set up initial conditions for gas cells */
  /* --------------------------------------- */
  double hydrogen_fraction_by_mass = 1.00;
  double gas_density = const_mh; /* corresponds to number density 1 cm^-3 */
  double T = 100.; /* K */

  /* Initial conditions for radiation */
  /* See README for details */
  double T_blackbody = 1e5; /* K */
  double frequency_bins[4] = {0., 3.288e15, 5.945e15, 13.157e15}; /* Hz */
  double fixed_luminosity[4] = {0., 1.350e+01, 2.779e+01, 6.152e+00}; /* erg / cm^2 / s */
  double fixed_radiation_density_field[4] = {0., 0., 0., 0.};
  for (int g = 0; g < RT_NGROUPS; g++){
    fixed_radiation_density_field[g] = fixed_luminosity[g] / const_speed_light_c;
  }

  /* Derived quantities from ICs */
  /* --------------------------- */

  /* Assuming fully neutral gas */
  double mu_init = 1.;

  double internal_energy = T * const_kboltz / (const_mh * mu_init * (const_adiabatic_index - 1.));

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

  /* estimate the mean weight based on the densities */
  double mean_weight = mean_weight_from_densities(
      gas_density, nHI, nHII, nHeI, nHeII, nHeIII, ne, mass_units);

  gr_float HI_density = nHI * (const_mh / mass_units);
  gr_float HII_density = nHII * (const_mh / mass_units);
  gr_float HeI_density = nHeI * (4 * const_mh / mass_units);
  gr_float HeII_density = nHeII * (4 * const_mh / mass_units);
  gr_float HeIII_density = nHeIII * (4 * const_mh / mass_units);
  /* !! this is the convention adopted by Grackle
   * !! e_density is the electron density * nh/ne */
  gr_float e_density = ne * (const_mh / mass_units);

  /* Store them all in a single array for simplicity. */
  gr_float species_densities[12] = {
      HI_density, HII_density, HeI_density, HeII_density, HeIII_density,
      e_density,  0.,          0.,          0.,           0.,
      0.,         0.};

  /* Get photon cross sections and mean energies */
  /* ------------------------------------------- */
  double **cse = malloc(RT_NGROUPS * sizeof(double *));
  double **csn = malloc(RT_NGROUPS * sizeof(double *));
  double mean_photon_energies[RT_NGROUPS];
  for (int group = 0; group < RT_NGROUPS; group++) {
    cse[group] = malloc(RT_NIONIZING_SPECIES * sizeof(double));
    csn[group] = malloc(RT_NIONIZING_SPECIES * sizeof(double));
    mean_photon_energies[group] = 0.;
  }

  get_cross_sections(T_blackbody, frequency_bins, cse, csn, mean_photon_energies);

  for (int g = 0; g < RT_NGROUPS; g++){
    printf("group %d cse=%.3e csn=%.3e mean=%.3e\n", g, cse[g][0], csn[g][0], mean_photon_energies[g]);
  }


  /* Grackle behaviour setup */
  /* ----------------------- */

  int UVbackground = 0;         /* toogle on/off the UV background */
  int primordial_chemistry = 1; /* choose the chemical network */
  int use_radiative_cooling = 1;
  int use_radiative_transfer = 0;
  char *grackle_data_file = "";

  /* Store them all in a single array for simplicity. */
  gr_float interaction_rates[5] = {0., 0., 0., 0., 0};

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
                          UVbackground, grackle_data_file, use_radiative_cooling, 
                          use_radiative_transfer, hydrogen_fraction_by_mass);

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

  /* Additional arrays to store temperature and mean molecular weights
   * of each cell. */
  gr_float *temperature = malloc(FIELD_SIZE * sizeof(gr_float));
  gr_float *mu = malloc(FIELD_SIZE * sizeof(gr_float));

  /* Write headers */
  /* ------------- */

  if (verbose) {
    printf("%15s%15s%15s%15s%15s%15s%15s%15s\n",
           "Initial setup: ", "Temperature", "nHI", "nHII", "nHeI", "nHeII",
           "nHeIII", "ne");
    printf("%15s%15g%15g%15g%15g%15g%15g%15g\n\n", "Initial setup: ", T, nHI,
           nHII, nHeI, nHeII, nHeIII, ne);
  }

  printf("%7s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n", "step", "Time [yr]",
         "dt [yr]", "Temperature [K]", "Mean Mol. W. [1]", "Tot dens. [IU]",
         "HI dens. [IU]", "HII dens. [IU]", "HeI dens. [IU]", "HeII dens. [IU]",
         "HeIII dens. [IU]", "e- num. dens. [IU]");
  /* TODO: dt is wrong here */
  printf(
      "%6d %15.3e%15.3e%15.1f%15.3f%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e \n", 0, 
      t / const_yr * time_units, dt_max / const_yr * time_units, T, mean_weight,
      grackle_fields.density[0], grackle_fields.HI_density[0],
      grackle_fields.HII_density[0], grackle_fields.HeI_density[0],
      grackle_fields.HeII_density[0], grackle_fields.HeIII_density[0],
      grackle_fields.e_density[0]);

  /* write down what ICs you used into file */
  /* TODO: dt is wrong here */
  write_my_setup(fd, grackle_fields, grackle_chemistry_data, mass_units, length_units, velocity_units, dt_max, hydrogen_fraction_by_mass, gas_density, internal_energy);
  fprintf(fd, "#%14s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n", "Time [yr]",
          "dt [yr]", "Temperature [K]", "Mean Mol. W. [1]", "Tot dens. [IU]",
          "HI dens. [IU]", "HII dens. [IU]", "HeI dens. [IU]",
          "HeII dens. [IU]", "HeIII dens. [IU]", "e- num. dens. [IU]");
  /* TODO: dt is wrong here */
  fprintf(
      fd,
      "%15.3e%15.3e%15.1f%15.3f%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e \n",
      t / const_yr * time_units, dt_max / const_yr * time_units, T, mean_weight,
      grackle_fields.density[0], grackle_fields.HI_density[0],
      grackle_fields.HII_density[0], grackle_fields.HeI_density[0],
      grackle_fields.HeII_density[0], grackle_fields.HeIII_density[0],
      grackle_fields.e_density[0]);

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/

  int step = 0;
  while (t < tend) {


    /* Set up radiation fields, and compute the resulting interaction
     * rates depending on the simulation time. */
    gr_float interaction_rates[5] = {0., 0., 0., 0., 0.};
    if (t / const_yr * time_units < 8.5e6) {
      /* below 0.5 Myr, we heat. */
      /* double test = ; */
      /* RT_HI_ionization_rate = 1e-22; */
      /* RT_HeI_ionization_rate = 0.; */
      /* RT_HeII_ionization_rate = 0.; */
      /* RT_H2_dissociation_rate = 0.; */
      /* RT_heating_rate = 1e-22; */

      /* RT_HI_ionization_rate = 1.6e-06; [> in s^-1 <] */
      /* RT_HeI_ionization_rate = 0; */
      /* RT_HeII_ionization_rate = 0; */
      /* RT_H2_dissociation_rate = 0; */
      /* RT_heating_rate = 5.191e-17 * 0.3; [> in erg/s/cm3 <] */

      /* printf("heating\n"); */

      double radiation_energy_density[RT_NGROUPS];
      radiation_energy_density[0] = fixed_radiation_density_field[0];
      radiation_energy_density[1] = fixed_radiation_density_field[1];
      radiation_energy_density[2] = fixed_radiation_density_field[2];
      radiation_energy_density[3] = fixed_radiation_density_field[3];

      gr_float species_densities[6];
      species_densities[0] = grackle_fields.HI_density[0];
      species_densities[1] = grackle_fields.HII_density[0];
      species_densities[2] = grackle_fields.HeI_density[0];
      species_densities[3] = grackle_fields.HeII_density[0];
      species_densities[4] = grackle_fields.HeIII_density[0];
      species_densities[5] = grackle_fields.e_density[0];

      get_interaction_rates(radiation_energy_density, 
                            species_densities,
                            cse, csn, mean_photon_energies,
                            interaction_rates);
    }

    for (int i = 0; i < FIELD_SIZE; i++){
      printf("%.4e\n", interaction_rates[0]);
      grackle_fields.RT_heating_rate[i] = interaction_rates[0];
      grackle_fields.RT_HI_ionization_rate[i] = interaction_rates[1];
      grackle_fields.RT_HeI_ionization_rate[i] = interaction_rates[2];
      grackle_fields.RT_HeII_ionization_rate[i] = interaction_rates[3];
      grackle_fields.RT_H2_dissociation_rate[i] = interaction_rates[4];
    }

    /* gr_float cooling_time[FIELD_SIZE]; */
    /* if (calculate_cooling_time(&grackle_units_data, &grackle_fields, */
    /*                            cooling_time) == 0) { */
    /*   fprintf(stderr, "Error in calculate_cooling_time.\n"); */
    /*   return 0; */
    /* } */
    /* gr_float dt = dt_max; */
    /* for (int i = 0; i < FIELD_SIZE; i++){ */
    /*   [> If we're heating, the cooling time will be negative <] */
    /*   if (fabs(cooling_time[i]) < dt) dt = fabs(cooling_time[i]); */
    /* } */
    /*  */
    gr_float dt = dt_max;
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

    if (mean_weight_local_like_grackle(&grackle_chemistry_data, &grackle_fields,
                                       mu) != SUCCESS) {
      fprintf(stderr, "Error in local_calculate_mean_weight.\n");
      return EXIT_FAILURE;
    }

    /* Calculate temperature. */
    if (calculate_temperature(&grackle_units_data, &grackle_fields,
                              temperature) == 0) {
      /* if (local_calculate_temperature(&grackle_chemistry_data,
       * grackle_chemistry_rates, &grackle_units_data, &grackle_fields,
       * temperature) == 0) { */
      fprintf(stderr, "Error in calculate_temperature.\n");
      return EXIT_FAILURE;
    }

    printf(
        "%6d %15.3e%15.3e%15.1f%15.3f%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e \n",
        step, t / const_yr * time_units, dt / const_yr * time_units, temperature[0],
        mu[0], grackle_fields.density[0], grackle_fields.HI_density[0],
        grackle_fields.HII_density[0], grackle_fields.HeI_density[0],
        grackle_fields.HeII_density[0], grackle_fields.HeIII_density[0],
        grackle_fields.e_density[0]);

    if (step % output_frequency == 0)
      fprintf(fd,
              "%15.3e%15.3e%15.1f%15.3f%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e%15."
              "3e \n",
              t / const_yr * time_units, dt / const_yr * time_units,
              temperature[0], mu[0], grackle_fields.density[0],
              grackle_fields.HI_density[0], grackle_fields.HII_density[0],
              grackle_fields.HeI_density[0], grackle_fields.HeII_density[0],
              grackle_fields.HeIII_density[0], grackle_fields.e_density[0]);
  }

  fclose(fd);

  clean_up_fields(&grackle_fields);
  free(mu);
  free(temperature);
  free(cse);
  free(csn);

  return EXIT_SUCCESS;
}
