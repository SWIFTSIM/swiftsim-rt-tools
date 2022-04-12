/* ---------------------------------------------
 * In this example, we start with high internal
 * energies and a fully ionized gas, and just
 * let it cool without any RT.
 * --------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <grackle.h>

#include "constants.h"
#include "ionization_equilibrium.h"
#include "mean_molecular_weight.h"

#define FIELD_SIZE 1
#define GRIDDIM 1

/*********************************************************************
/ Grackle Setup Functions
*********************************************************************/

/* Set up the units for grackle. */
void setup_grackle_units(code_units *grackle_units_data, double density_units,
                         double length_units, double time_units) {

  grackle_units_data->comoving_coordinates = 0; /* no cosmo */
  grackle_units_data->density_units = density_units;
  grackle_units_data->length_units = length_units;
  grackle_units_data->time_units = time_units;
  grackle_units_data->velocity_units =
      grackle_units_data->length_units / grackle_units_data->time_units;
  grackle_units_data->a_units = 1.0; /* units for the expansion factor */
  /* Set expansion factor to 1 for non-cosmological simulation-> */
  grackle_units_data->a_value = 1.;

  /* set temperature units */
  /* double temperature_units = */
  /*     const_mh / const_kboltz * */
  /*     pow(grackle_units_data.a_units * grackle_units_data.length_units / */
  /*             grackle_units_data.time_units, */
  /* 2); */
}

/* Set parameters for the grackle chemistry object */
void setup_grackle_chemistry(chemistry_data *grackle_chemistry_data,
                             int primordial_chemistry, int UVbackground,
                             char *grackle_data_file,
                             int use_radiative_transfer,
                             float hydrogen_fraction_by_mass) {
  /* Set parameter values for chemistry. */
  grackle_chemistry_data->use_grackle = 1;            /* chemistry on */
  grackle_chemistry_data->with_radiative_cooling = 1; /* cooling on */
  grackle_chemistry_data->primordial_chemistry = primordial_chemistry;
  grackle_chemistry_data->dust_chemistry = 0;
  grackle_chemistry_data->metal_cooling = 0;
  grackle_chemistry_data->UVbackground = UVbackground;
  grackle_chemistry_data->grackle_data_file = grackle_data_file;
  grackle_chemistry_data->use_radiative_transfer = use_radiative_transfer;
  grackle_chemistry_data->HydrogenFractionByMass = hydrogen_fraction_by_mass;
  grackle_chemistry_data->Gamma =
      const_adiabatic_index; /* defined in const->h */
}

/* Set up grackle gas fields. */
void setup_grackle_fields(grackle_field_data *grackle_fields,
                          gr_float species_densities[12],
                          gr_float interaction_rates[5], double gas_density,
                          double internal_energy) {

  /* Set grid dimension and size.
   * grid_start and grid_end are used to ignore ghost zones. */
  grackle_fields->grid_rank = GRIDDIM;
  grackle_fields->grid_dimension =
      malloc(grackle_fields->grid_rank * sizeof(int));
  grackle_fields->grid_start = malloc(grackle_fields->grid_rank * sizeof(int));
  grackle_fields->grid_end = malloc(grackle_fields->grid_rank * sizeof(int));
  /* used only for H2 self-shielding approximation */
  /* grackle_fields->grid_dx = 0.0; */

  /* for (int i = 0; i < GRIDDIM; i++) { */
  /*   [> the active dimension not including ghost zones. <] */
  /*   grackle_fields->grid_dimension[i] = FIELD_SIZE;  */
  /*   grackle_fields->grid_start[i] = 0; */
  /*   grackle_fields->grid_end[i] = 0; */
  /* } */
  grackle_fields->grid_dimension[0] = FIELD_SIZE;
  grackle_fields->grid_start[0] = 0;
  grackle_fields->grid_end[0] = FIELD_SIZE - 1;

  /* Set initial quantities */
  grackle_fields->density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->internal_energy = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->x_velocity = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->y_velocity = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->z_velocity = malloc(FIELD_SIZE * sizeof(gr_float));
  /* for primordial_chemistry >= 1 */
  grackle_fields->HI_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->HII_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->HeI_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->HeII_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->HeIII_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->e_density = malloc(FIELD_SIZE * sizeof(gr_float));
  /* for primordial_chemistry >= 2 */
  grackle_fields->HM_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->H2I_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->H2II_density = malloc(FIELD_SIZE * sizeof(gr_float));
  /* for primordial_chemistry >= 3 */
  grackle_fields->DI_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->DII_density = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->HDI_density = malloc(FIELD_SIZE * sizeof(gr_float));
  /* for metal_cooling = 1 */
  grackle_fields->metal_density = malloc(FIELD_SIZE * sizeof(gr_float));

  /* volumetric heating rate (provide in units [erg s^-1 cm^-3]) */
  grackle_fields->volumetric_heating_rate =
      malloc(FIELD_SIZE * sizeof(gr_float));
  /* specific heating rate (provide in units [egs s^-1 g^-1] */
  grackle_fields->specific_heating_rate = malloc(FIELD_SIZE * sizeof(gr_float));

  /* radiative transfer ionization / dissociation rate fields (provide in units
   * [1/s]) */
  grackle_fields->RT_HI_ionization_rate = malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->RT_HeI_ionization_rate =
      malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->RT_HeII_ionization_rate =
      malloc(FIELD_SIZE * sizeof(gr_float));
  grackle_fields->RT_H2_dissociation_rate =
      malloc(FIELD_SIZE * sizeof(gr_float));
  /* radiative transfer heating rate field (provide in units [erg s^-1 cm^-3])
   */
  grackle_fields->RT_heating_rate = malloc(FIELD_SIZE * sizeof(gr_float));

  for (int i = 0; i < FIELD_SIZE; i++) {

    /* initial density */
    grackle_fields->density[i] = gas_density;

    /* initial internal energy using initial temperature */
    grackle_fields->internal_energy[i] = internal_energy;
    /* T / (grackle_chemistry_data.Gamma - 1.0) / mean_weight /
     * temperature_units; */

    grackle_fields->HI_density[i] = species_densities[0];
    grackle_fields->HII_density[i] = species_densities[1];
    grackle_fields->HeI_density[i] = species_densities[2];
    grackle_fields->HeII_density[i] = species_densities[3];
    grackle_fields->HeIII_density[i] = species_densities[4];
    grackle_fields->e_density[i] =
        species_densities[5]; /* electron density*mh/me */

    grackle_fields->HM_density[i] = species_densities[6];
    grackle_fields->H2I_density[i] = species_densities[7];
    grackle_fields->H2II_density[i] = species_densities[8];
    grackle_fields->DI_density[i] = species_densities[9];
    grackle_fields->DII_density[i] = species_densities[10];
    grackle_fields->HDI_density[i] = species_densities[11];

    /* solar metallicity */
    grackle_fields->metal_density[i] = 0.0;
    /* grackle_chemistry_data.SolarMetalFractionByMass *
     * grackle_fields->density[i]; */

    grackle_fields->x_velocity[i] = 0.0;
    grackle_fields->y_velocity[i] = 0.0;
    grackle_fields->z_velocity[i] = 0.0;

    grackle_fields->volumetric_heating_rate[i] = 0.0;
    grackle_fields->specific_heating_rate[i] = 0.0;

    grackle_fields->RT_HI_ionization_rate[i] = interaction_rates[0];
    grackle_fields->RT_HeI_ionization_rate[i] = interaction_rates[1];
    grackle_fields->RT_HeII_ionization_rate[i] = interaction_rates[2];
    grackle_fields->RT_H2_dissociation_rate[i] = interaction_rates[3];
    grackle_fields->RT_heating_rate[i] = interaction_rates[4];
  }
}

/* Deallocate fields when you're done. */
void clean_up_fields(grackle_field_data *grackle_fields) {

  free(grackle_fields->grid_dimension);

  free(grackle_fields->grid_start);
  free(grackle_fields->grid_end);
  free(grackle_fields->density);
  free(grackle_fields->internal_energy);
  free(grackle_fields->x_velocity);
  free(grackle_fields->y_velocity);
  free(grackle_fields->z_velocity);
  free(grackle_fields->HI_density);
  free(grackle_fields->HII_density);
  free(grackle_fields->HeI_density);
  free(grackle_fields->HeII_density);
  free(grackle_fields->HeIII_density);
  free(grackle_fields->e_density);
  free(grackle_fields->HM_density);
  free(grackle_fields->H2I_density);
  free(grackle_fields->H2II_density);
  free(grackle_fields->DI_density);
  free(grackle_fields->DII_density);
  free(grackle_fields->HDI_density);
  free(grackle_fields->metal_density);
  free(grackle_fields->volumetric_heating_rate);
  free(grackle_fields->specific_heating_rate);

  free(grackle_fields->RT_HI_ionization_rate);
  free(grackle_fields->RT_HeI_ionization_rate);
  free(grackle_fields->RT_HeII_ionization_rate);
  free(grackle_fields->RT_H2_dissociation_rate);
  free(grackle_fields->RT_heating_rate);
}

/*********************************************************************
/ Main
*********************************************************************/

int main() {

  /******************************************************************
   * Set up initial conditions and runtime parameters.
   *****************************************************************/

  /* Print some extra data to screen? */
  int verbose = 1;
  /* output file */
  FILE *fd = fopen("out.dat", "w");
  /* output frequency  in number of steps */
  int output_frequency = 8;

  /* Define units : use the same as internal units for swift */
  /* ------------------------------------------------------- */
  double mass_units = 1.99848e43;
  double length_units = 3.08567758e21;
  double velocity_units = 1e5;

  double density_units =
      mass_units / length_units / length_units / length_units;
  double time_units = length_units / velocity_units;

  /* Time integration variables */
  /* -------------------------- */
  double t = 0.;
  double dt = 4.882813e-05; /* in internal units. Copy this from swift output */
  double tinit = 1e3;       /* in yr; will be converted later */
  double tend = 1e8;        /* in yr; will be converted later */
  t = tinit * const_yr / time_units;   /* yr to code units */
  tend = tend * const_yr / time_units; /* yr to code units */

  /* Set up initial conditions for gas cells */
  /* --------------------------------------- */
  double hydrogen_fraction_by_mass = 0.76;
  /* Use solution of swift's output. This is in internal units already. */
  double gas_density = 0.00024633363;
  double internal_energy = 21201.879;

  /* Derived quantities from ICs */
  /* --------------------------- */

  /* Assuming fully ionized gas */
  double mu_init = mean_molecular_weight_from_mass_fractions(
      0., hydrogen_fraction_by_mass, 0., 0., (1. - hydrogen_fraction_by_mass));
  double internal_energy_cgs =
      internal_energy * length_units * length_units / time_units / time_units;

  double T = internal_energy_cgs * (const_adiabatic_index - 1) * mu_init *
             const_mh / const_kboltz;
  if (verbose)
    printf("Initial setup: u_cgs %g T_cgs %g\n", internal_energy_cgs, T);

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

  /* Grackle behaviour setup */
  /* ----------------------- */

  int UVbackground = 0;         /* toogle on/off the UV background */
  int primordial_chemistry = 1; /* choose the chemical network */

  int use_radiative_transfer = 0;
  gr_float RT_HI_ionization_rate = 0.;
  gr_float RT_HeI_ionization_rate = 0.;
  gr_float RT_HeII_ionization_rate = 0.;
  gr_float RT_H2_dissociation_rate = 0.;
  gr_float RT_heating_rate = 0.;
  char *grackle_data_file = "";

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
    return EXIT_FAILURE;
  }

  /* Set parameter values for chemistry. */
  setup_grackle_chemistry(&grackle_chemistry_data, primordial_chemistry,
                          UVbackground, grackle_data_file,
                          use_radiative_transfer, hydrogen_fraction_by_mass);

  /* Initialize the chemistry_data_storage object to be able to use local
   * functions */
  /* chemistry_data_storage *grackle_chemistry_rates; */
  /* grackle_chemistry_rates = malloc(sizeof(chemistry_data_storage)); */
  /* if (_initialize_chemistry_data(&grackle_chemistry_data,
   * grackle_chemistry_rates, &grackle_units_data) == 0) { */
  /*   fprintf(stderr, "Error in initialize_chemistry_data.\n"); */
  /*   return EXIT_FAILURE; */
  /* } */
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

  if (verbose) {
    printf("%15s%15s%15s%15s%15s%15s%15s%15s\n",
           "Initial setup: ", "Temperature", "nHI", "nHII", "nHeI", "nHeII",
           "nHeIII", "ne");
    printf("%15s%15g%15g%15g%15g%15g%15g%15g\n\n", "Initial setup: ", T, nHI,
           nHII, nHeI, nHeII, nHeIII, ne);
  }

  printf("%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n", "Time [yr]",
         "dt [yr]", "Temperature [K]", "Mean Mol. W. [1]", "Tot dens. [IU]",
         "HI dens. [IU]", "HII dens. [IU]", "HeI dens. [IU]", "HeII dens. [IU]",
         "HeIII dens. [IU]", "e- num. dens. [IU]");
  printf(
      "%15.3e%15.3e%15.1f%15.3f%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e \n",
      t / const_yr * time_units, dt / const_yr * time_units, T, mean_weight,
      grackle_fields.density[0], grackle_fields.HI_density[0],
      grackle_fields.HII_density[0], grackle_fields.HeI_density[0],
      grackle_fields.HeII_density[0], grackle_fields.HeIII_density[0],
      grackle_fields.e_density[0]);

  /* write down what ICs you used into file */
  fprintf(fd, "# Result file created using grackle standalone program.\n");
  fprintf(fd, "# mass units used: %.6g [g]\n", mass_units);
  fprintf(fd, "# length units used: %.6g [cm]\n", length_units);
  fprintf(fd, "# velocity units units used: %.6g [cm/s]\n", velocity_units);
  fprintf(fd, "# dt used: %.6g [internal units]\n", dt);
  fprintf(fd, "# hydrogen mass fraction used: %.6g\n",
          hydrogen_fraction_by_mass);
  fprintf(fd, "# gas density used: %.6g [internal units]\n", gas_density);
  fprintf(fd, "# inital internal density used: %.6g [internal units]\n",
          internal_energy);
  fprintf(fd, "# Grackle parameters:\n");
  fprintf(fd, "# grackle_chemistry_data.use_grackle = %d\n",
          grackle_chemistry_data.use_grackle);
  fprintf(fd, "# grackle_chemistry_data.with_radiative_cooling %d\n",
          grackle_chemistry_data.with_radiative_cooling);
  fprintf(fd, "# grackle_chemistry_data.primordial_chemistry = %d\n",
          grackle_chemistry_data.primordial_chemistry);
  fprintf(fd, "# grackle_chemistry_data.dust_chemistry = %d\n",
          grackle_chemistry_data.dust_chemistry);
  fprintf(fd, "# grackle_chemistry_data.metal_cooling = %d\n",
          grackle_chemistry_data.metal_cooling);
  fprintf(fd, "# grackle_chemistry_data.UVbackground = %d\n",
          grackle_chemistry_data.UVbackground);
  fprintf(fd, "# grackle_chemistry_data.grackle_data_file = %s\n",
          grackle_chemistry_data.grackle_data_file);
  fprintf(fd, "# grackle_chemistry_data.use_radiative_transfer = %d\n",
          grackle_chemistry_data.use_radiative_transfer);
  fprintf(fd, "# grackle_chemistry_data.HydrogenFractionByMass = %.3g\n",
          grackle_chemistry_data.HydrogenFractionByMass);
  fprintf(fd, "# grackle_chemistry_data.Gamma = %.6g\n",
          grackle_chemistry_data.Gamma);
  fprintf(fd, "# Grackle field data:\n");

  fprintf(fd, "# grackle_fields.density = %.6g\n", grackle_fields.density[0]);

  fprintf(fd, "# grackle_fields.internal_energy = %.6g\n",
          grackle_fields.internal_energy[0]);

  fprintf(fd, "# grackle_fields.HI_density = %.6g\n",
          grackle_fields.HI_density[0]);

  fprintf(fd, "# grackle_fields.HII_density = %.6g\n",
          grackle_fields.HII_density[0]);
  fprintf(fd, "# grackle_fields.HeI_density = %.6g\n",
          grackle_fields.HeI_density[0]);
  fprintf(fd, "# grackle_fields.HeII_density = %.6g\n",
          grackle_fields.HeII_density[0]);
  fprintf(fd, "# grackle_fields.HeIII_density = %.6g\n",
          grackle_fields.HeIII_density[0]);
  fprintf(fd, "# grackle_fields.e_density = %.6g\n",
          grackle_fields.e_density[0]);

  fprintf(fd, "# grackle_fields.HM_density  = %.6g\n",
          grackle_fields.HM_density[0]);
  fprintf(fd, "# grackle_fields.H2I_density = %.6g\n",
          grackle_fields.H2I_density[0]);
  fprintf(fd, "# grackle_fields.H2II_density = %.6g\n",
          grackle_fields.H2II_density[0]);
  fprintf(fd, "# grackle_fields.DI_density = %.6g\n",
          grackle_fields.DI_density[0]);
  fprintf(fd, "# grackle_fields.DII_density = %.6g\n",
          grackle_fields.DII_density[0]);
  fprintf(fd, "# grackle_fields.HDI_density = %.6g\n",
          grackle_fields.HDI_density[0]);

  fprintf(fd, "# grackle_fields.metal_density = %.6g\n",
          grackle_fields.metal_density[0]);

  fprintf(fd, "# grackle_fields.x_velocity = %.6g\n",
          grackle_fields.x_velocity[0]);
  fprintf(fd, "# grackle_fields.y_velocity = %.6g\n",
          grackle_fields.y_velocity[0]);
  fprintf(fd, "# grackle_fields.z_velocity = %.6g\n",
          grackle_fields.z_velocity[0]);

  fprintf(fd, "# grackle_fields.volumetric_heating_rate = %.6g\n",
          grackle_fields.volumetric_heating_rate[0]);
  fprintf(fd, "# grackle_fields.specific_heating_rate = %.6g\n",
          grackle_fields.specific_heating_rate[0]);

  fprintf(fd, "# grackle_fields.RT_HI_ionization_rate = %.6g\n",
          grackle_fields.RT_HI_ionization_rate[0]);
  fprintf(fd, "# grackle_fields.RT_HeI_ionization_rate = %.6g\n",
          grackle_fields.RT_HeI_ionization_rate[0]);
  fprintf(fd, "# grackle_fields.RT_HeII_ionization_rate = %.6g\n",
          grackle_fields.RT_HeII_ionization_rate[0]);
  fprintf(fd, "# grackle_fields.RT_H2_dissociation_rate = %.6g\n",
          grackle_fields.RT_H2_dissociation_rate[0]);
  fprintf(fd, "# grackle_fields.RT_heating_rate = %.6g\n",
          grackle_fields.RT_heating_rate[0]);

  fprintf(fd, "#%14s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n", "Time [yr]",
          "dt [yr]", "Temperature [K]", "Mean Mol. W. [1]", "Tot dens. [IU]",
          "HI dens. [IU]", "HII dens. [IU]", "HeI dens. [IU]",
          "HeII dens. [IU]", "HeIII dens. [IU]", "e- num. dens. [IU]");
  fprintf(
      fd,
      "%15.3e%15.3e%15.1f%15.3f%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e \n",
      t / const_yr * time_units, dt / const_yr * time_units, T, mean_weight,
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
        "%15.3e%15.3e%15.1f%15.3f%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e%15.3e \n",
        t / const_yr * time_units, dt / const_yr * time_units, temperature[0],
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

  return EXIT_SUCCESS;
}
