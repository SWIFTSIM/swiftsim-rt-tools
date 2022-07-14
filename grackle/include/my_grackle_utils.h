#ifndef MY_GRACKLE_UTILS_H
#define MY_GRACKLE_UTILS_H

/*-----------------------------------------------------
 * Contains generic grackle utilities, e.g. for setup,
 * printouts, etc.
 *--------------------------------------------------- */

#include <grackle.h>
#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "mean_molecular_weight.h"

/**
 * @brief Set up the units for grackle.
 **/
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

/**
 * @brief Set parameters for the grackle chemistry object
 **/
void setup_grackle_chemistry(chemistry_data *grackle_chemistry_data,
                             int primordial_chemistry, int UVbackground,
                             char *grackle_data_file, int use_radiative_cooling,
                             int use_radiative_transfer,
                             float hydrogen_fraction_by_mass) {
  /* Set parameter values for chemistry. */
  grackle_chemistry_data->use_grackle = 1; /* chemistry on */
  grackle_chemistry_data->with_radiative_cooling = use_radiative_cooling;
  grackle_chemistry_data->primordial_chemistry = primordial_chemistry;
  grackle_chemistry_data->dust_chemistry = 0;
  grackle_chemistry_data->metal_cooling = 0;
  grackle_chemistry_data->UVbackground = UVbackground;
  grackle_chemistry_data->grackle_data_file = grackle_data_file;
  grackle_chemistry_data->use_radiative_transfer = use_radiative_transfer;
  grackle_chemistry_data->HydrogenFractionByMass = hydrogen_fraction_by_mass;
  grackle_chemistry_data->Gamma =
      const_adiabatic_index; /* defined in const.h */
}

/**
 * @brief Set up grackle gas fields.
 **/
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

    grackle_fields->RT_heating_rate[i] = interaction_rates[0];
    grackle_fields->RT_HI_ionization_rate[i] = interaction_rates[1];
    grackle_fields->RT_HeI_ionization_rate[i] = interaction_rates[2];
    grackle_fields->RT_HeII_ionization_rate[i] = interaction_rates[3];
    grackle_fields->RT_H2_dissociation_rate[i] = interaction_rates[4];
  }
}

/**
 * @brief Deallocate fields when you're done.
 **/
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

/**
 * @brief Dump the setup used to generate the example
 **/
void write_my_setup(FILE *fd, grackle_field_data grackle_fields,
                    chemistry_data grackle_chemistry_data, double mass_units,
                    double length_units, double velocity_units, double dt,
                    double hydrogen_fraction_by_mass, double gas_density,
                    double internal_energy) {
  fprintf(fd, "# Result file created using grackle standalone program.\n");
  fprintf(fd, "# mass units used: %.6g [g]\n", mass_units);
  fprintf(fd, "# length units used: %.6g [cm]\n", length_units);
  fprintf(fd, "# velocity units units used: %.6g [cm/s]\n", velocity_units);
  fprintf(fd, "# dt used: %.6g [internal units]\n", dt);
  fprintf(fd, "# hydrogen mass fraction used: %.6g\n",
          hydrogen_fraction_by_mass);
  fprintf(fd, "# gas density used: %.6g [internal units]\n", gas_density);
  fprintf(fd, "# inital internal energy used: %.6g [internal units]\n",
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
}

/**
 * @brief write header to a file/stdout
 **/
void write_header(FILE *fd) {

  fprintf(fd, "#%8s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n",
          "step", "Time [yr]", "dt [yr]", "Temperature [K]", "Mean M Wgt [1]",
          "Tot dens [IU]", "HI dens [IU]", "HII dens [IU]", "HeI dens [IU]",
          "HeII dens [IU]", "HeIII dens [IU]", "e- n. dens [IU]");
}
/**
 * @brief write the current state of a field with index i to a file/stdout
 **/
void write_timestep(FILE *fd, grackle_field_data *grackle_fields,
                    code_units *grackle_units_data,
                    chemistry_data *grackle_chemistry_data, int field_index,
                    double t, double dt, double time_units, int step) {

  /* Additional arrays to store temperature and mean molecular weights
   * of each cell. */
  gr_float temperature[FIELD_SIZE];
  gr_float mu[FIELD_SIZE];
  for (int i = 0; i < FIELD_SIZE; i++) {
    temperature[i] = 0;
    mu[i] = 0;
  }

  /* Grab temperature and mean molecular weights. */
  /* Calculate temperature. */
  if (calculate_temperature(grackle_units_data, grackle_fields, temperature) ==
      0) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    abort();
  }

  if (mean_weight_local_like_grackle(grackle_chemistry_data, grackle_fields,
                                     mu) != SUCCESS) {
    fprintf(stderr, "Error in local_calculate_mean_weight.\n");
    abort();
  }

  fprintf(fd,
          "%9d %15.3e %15.3e %15.3e %15.3e %15.3e %15.3e %15.3e %15.3e %15.3e "
          "%15.3e %15.3e\n",
          step, t / const_yr * time_units, dt / const_yr * time_units,
          temperature[field_index], mu[field_index],
          grackle_fields->density[field_index],
          grackle_fields->HI_density[field_index],
          grackle_fields->HII_density[field_index],
          grackle_fields->HeI_density[field_index],
          grackle_fields->HeII_density[field_index],
          grackle_fields->HeIII_density[field_index],
          grackle_fields->e_density[field_index]);
}
#endif
