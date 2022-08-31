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
  grackle_units_data->a_units = 1.0;
  grackle_units_data->a_value = 1.;

  /* Set velocity units */
  set_velocity_units(grackle_units_data);
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

  int *dimension = malloc(3 * sizeof(int));
  int *start = malloc(3 * sizeof(int));
  int *end = malloc(3 * sizeof(int));

  /* Set grid dimension and size.
   * grid_start and grid_end are used to ignore ghost zones. */
  grackle_fields->grid_rank = 3;
  grackle_fields->grid_dimension = dimension;
  grackle_fields->grid_start = start;
  grackle_fields->grid_end = end;
  /* used only for H2 self-shielding approximation */
  /* grackle_fields->grid_dx = 0.0; */

  /* NOTE: if you're trying to simplify this, you MUST allocate GRIDDIM = 3
   * and grid_dimension, grid_start, grid_end with at least 3D as well.
   * Otherwise, grackle will cause segfaults because they do pointer arithmetics
   * assuming 3 dimensions internally. */
  for (int i = 0; i < 3; i++) {
    /* the active dimension not including ghost zones. */
    grackle_fields->grid_dimension[i] = 1;
    grackle_fields->grid_start[i] = 0;
    grackle_fields->grid_end[i] = 0;
  }
  grackle_fields->grid_dimension[0] = FIELD_SIZE;
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
 *
 * @param fp FILE pointer to write into.
 * @param grackle_fields grackle field data
 * @param grackle_chemistry_data grackle chemistry data.
 * @param mass_units mass units that convert internal units to cgs
 * @param length_units length units that convert internal units to cgs
 * @param velocity_units velocity units that convert internal units to cgs
 * @param dt current time step
 * @param hydrogen_fraction_by_mass hydrogen fraction by mass used.
 * @param gas_density gas density used. In internal units.
 * @param internal_energy internal energy used. In internal units.
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
  fprintf(fd, "# grackle_chemistry_data.CaseBRecombination = %d\n",
          grackle_chemistry_data.CaseBRecombination);
  fprintf(fd, "# grackle_chemistry_data.grackle_data_file = %s\n",
          grackle_chemistry_data.grackle_data_file);
  fprintf(fd, "# grackle_chemistry_data.use_radiative_transfer = %d\n",
          grackle_chemistry_data.use_radiative_transfer);
  fprintf(fd, "# grackle_chemistry_data.HydrogenFractionByMass = %.3g\n",
          grackle_chemistry_data.HydrogenFractionByMass);
  fprintf(fd, "# grackle_chemistry_data.Gamma = %.6g\n",
          grackle_chemistry_data.Gamma);
  fprintf(fd, "# Grackle field data:\n");

#define write_grackle_field(v)                                                 \
  if (grackle_fields.v != NULL)                                                \
  fprintf(fd, "# grackle_fields." #v " = %g\n", grackle_fields.v[0])

  write_grackle_field(density);
  write_grackle_field(internal_energy);
  write_grackle_field(HI_density);
  write_grackle_field(HII_density);
  write_grackle_field(HeI_density);
  write_grackle_field(HeII_density);
  write_grackle_field(HeIII_density);
  write_grackle_field(e_density);
  write_grackle_field(HM_density);
  write_grackle_field(H2I_density);
  write_grackle_field(H2II_density);
  write_grackle_field(DI_density);
  write_grackle_field(DII_density);
  write_grackle_field(HDI_density);
  write_grackle_field(metal_density);
  write_grackle_field(x_velocity);
  write_grackle_field(y_velocity);
  write_grackle_field(z_velocity);
  write_grackle_field(volumetric_heating_rate);
  write_grackle_field(specific_heating_rate);
  write_grackle_field(RT_HI_ionization_rate);
  write_grackle_field(RT_HeI_ionization_rate);
  write_grackle_field(RT_HeII_ionization_rate);
  write_grackle_field(RT_H2_dissociation_rate);
  write_grackle_field(RT_heating_rate);
}

/**
 * @brief write header to a file/stdout
 **/
void write_header(FILE *fd) {

  fprintf(fd,
          "#%8s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n",
          "step", "Time [yr]", "dt [yr]", "Temperature [K]", "Mean M Wgt [1]",
          "Tot dens [IU]", "HI dens [IU]", "HII dens [IU]", "HeI dens [IU]",
          "HeII dens [IU]", "HeIII dens [IU]", "e- n. dens [IU]",
          "IntrnEnerg [IU]");
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
    temperature[i] = 0.;
    mu[i] = 0.;
  }

  /* Grab temperature and mean molecular weights. */
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
          "%15.3e "
          "%15.3e %15.3e\n",
          step, t / const_yr * time_units, dt / const_yr * time_units,
          temperature[field_index], mu[field_index],
          grackle_fields->density[field_index],
          grackle_fields->HI_density[field_index],
          grackle_fields->HII_density[field_index],
          grackle_fields->HeI_density[field_index],
          grackle_fields->HeII_density[field_index],
          grackle_fields->HeIII_density[field_index],
          grackle_fields->e_density[field_index],
          grackle_fields->internal_energy[field_index]);
}

/**
 * @brief Write out all available grackle field data for a given index
 * and setup to a file.
 *
 * @param fp FILE pointer to write into.
 * @param grackle_fields grackle field data
 * @param grackle_chemistry_data grackle chemistry data.
 * @param grackle_units units used by grackle
 * @param field_index grackle field index to print out.
 **/
void print_grackle_setup_and_field(FILE *fp, grackle_field_data grackle_fields,
                                   chemistry_data *grackle_chemistry_data,
                                   code_units grackle_units, int field_index) {

  fprintf(fp, "Grackle chemistry parameters:\n");

  fprintf(fp, "use_grackle                       = %d\n",
          grackle_chemistry_data->use_grackle);
  fprintf(fp, "with_radiative_cooling            = %d\n",
          grackle_chemistry_data->with_radiative_cooling);
  fprintf(fp, "primordial_chemistry              = %d\n",
          grackle_chemistry_data->primordial_chemistry);
  fprintf(fp, "dust_chemistry                    = %d\n",
          grackle_chemistry_data->dust_chemistry);
  fprintf(fp, "metal_cooling                     = %d\n",
          grackle_chemistry_data->metal_cooling);
  fprintf(fp, "UVbackground                      = %d\n",
          grackle_chemistry_data->UVbackground);
  fprintf(fp, "grackle_data_file                 = %s\n",
          grackle_chemistry_data->grackle_data_file);
  fprintf(fp, "cmb_temperature_floor             = %d\n",
          grackle_chemistry_data->cmb_temperature_floor);
  fprintf(fp, "Gamma                             = %g\n",
          grackle_chemistry_data->Gamma);
  fprintf(fp, "h2_on_dust                        = %d\n",
          grackle_chemistry_data->h2_on_dust);
  fprintf(fp, "use_dust_density_field            = %d\n",
          grackle_chemistry_data->use_dust_density_field);
  fprintf(fp, "dust_recombination_cooling        = %d\n",
          grackle_chemistry_data->dust_recombination_cooling);
  fprintf(fp, "photoelectric_heating             = %d\n",
          grackle_chemistry_data->photoelectric_heating);
  fprintf(fp, "photoelectric_heating_rate        = %g\n",
          grackle_chemistry_data->photoelectric_heating_rate);
  fprintf(fp, "use_isrf_field                    = %d\n",
          grackle_chemistry_data->use_isrf_field);
  fprintf(fp, "interstellar_radiation_field      = %g\n",
          grackle_chemistry_data->interstellar_radiation_field);
  fprintf(fp, "use_volumetric_heating_rate       = %d\n",
          grackle_chemistry_data->use_volumetric_heating_rate);
  fprintf(fp, "use_specific_heating_rate         = %d\n",
          grackle_chemistry_data->use_specific_heating_rate);
  fprintf(fp, "three_body_rate                   = %d\n",
          grackle_chemistry_data->three_body_rate);
  fprintf(fp, "cie_cooling                       = %d\n",
          grackle_chemistry_data->cie_cooling);
  fprintf(fp, "h2_optical_depth_approximation    = %d\n",
          grackle_chemistry_data->h2_optical_depth_approximation);
  fprintf(fp, "ih2co                             = %d\n",
          grackle_chemistry_data->ih2co);
  fprintf(fp, "ipiht                             = %d\n",
          grackle_chemistry_data->ipiht);
  fprintf(fp, "HydrogenFractionByMass            = %g\n",
          grackle_chemistry_data->HydrogenFractionByMass);
  fprintf(fp, "DeuteriumToHydrogenRatio          = %g\n",
          grackle_chemistry_data->DeuteriumToHydrogenRatio);
  fprintf(fp, "SolarMetalFractionByMass          = %g\n",
          grackle_chemistry_data->SolarMetalFractionByMass);
  fprintf(fp, "local_dust_to_gas_ratio           = %g\n",
          grackle_chemistry_data->local_dust_to_gas_ratio);
  fprintf(fp, "NumberOfTemperatureBins           = %d\n",
          grackle_chemistry_data->NumberOfTemperatureBins);
  fprintf(fp, "CaseBRecombination                = %d\n",
          grackle_chemistry_data->CaseBRecombination);
  fprintf(fp, "TemperatureStart                  = %g\n",
          grackle_chemistry_data->TemperatureStart);
  fprintf(fp, "TemperatureEnd                    = %g\n",
          grackle_chemistry_data->TemperatureEnd);
  fprintf(fp, "NumberOfDustTemperatureBins       = %d\n",
          grackle_chemistry_data->NumberOfDustTemperatureBins);
  fprintf(fp, "DustTemperatureStart              = %g\n",
          grackle_chemistry_data->DustTemperatureStart);
  fprintf(fp, "DustTemperatureEnd                = %g\n",
          grackle_chemistry_data->DustTemperatureEnd);
  fprintf(fp, "Compton_xray_heating              = %d\n",
          grackle_chemistry_data->Compton_xray_heating);
  fprintf(fp, "LWbackground_sawtooth_suppression = %d\n",
          grackle_chemistry_data->LWbackground_sawtooth_suppression);
  fprintf(fp, "LWbackground_intensity            = %g\n",
          grackle_chemistry_data->LWbackground_intensity);
  fprintf(fp, "UVbackground_redshift_on          = %g\n",
          grackle_chemistry_data->UVbackground_redshift_on);
  fprintf(fp, "UVbackground_redshift_off         = %g\n",
          grackle_chemistry_data->UVbackground_redshift_off);
  fprintf(fp, "UVbackground_redshift_fullon      = %g\n",
          grackle_chemistry_data->UVbackground_redshift_fullon);
  fprintf(fp, "UVbackground_redshift_drop        = %g\n",
          grackle_chemistry_data->UVbackground_redshift_drop);
  fprintf(fp, "cloudy_electron_fraction_factor   = %g\n",
          grackle_chemistry_data->cloudy_electron_fraction_factor);
  fprintf(fp, "use_radiative_transfer            = %d\n",
          grackle_chemistry_data->use_radiative_transfer);
  fprintf(fp, "radiative_transfer_coupled_rate_solver = %d\n",
          grackle_chemistry_data->radiative_transfer_coupled_rate_solver);
  fprintf(fp, "radiative_transfer_intermediate_step = %d\n",
          grackle_chemistry_data->radiative_transfer_intermediate_step);
  fprintf(fp, "radiative_transfer_hydrogen_only  = %d\n",
          grackle_chemistry_data->radiative_transfer_hydrogen_only);
  fprintf(fp, "self_shielding_method             = %d\n",
          grackle_chemistry_data->self_shielding_method);
  fprintf(fp, "H2_custom_shielding               = %d\n",
          grackle_chemistry_data->H2_custom_shielding);
  fprintf(fp, "H2_self_shielding                 = %d\n",
          grackle_chemistry_data->H2_self_shielding);

  fprintf(fp, "\nUnits:\n");
  fprintf(fp, "a_units               = %g\n", grackle_units.a_units);
  fprintf(fp, "a_value               = %g\n", grackle_units.a_value);
  fprintf(fp, "comoving_coordinates  = %d\n",
          grackle_units.comoving_coordinates);
  fprintf(fp, "density_units         = %g\n", grackle_units.density_units);
  fprintf(fp, "length_units          = %g\n", grackle_units.length_units);
  fprintf(fp, "time_units            = %g\n", grackle_units.time_units);
  fprintf(fp, "velocity_units        = %g\n", grackle_units.velocity_units);

#define rt_print_grackle_field(v)                                              \
  if (grackle_fields.v != NULL)                                                \
  fprintf(fp, "grackle_fields." #v " = %g\n", grackle_fields.v[field_index])

  fprintf(fp, "\nGrackle field data:\n");
  rt_print_grackle_field(density);
  rt_print_grackle_field(internal_energy);
  rt_print_grackle_field(HI_density);
  rt_print_grackle_field(HII_density);
  rt_print_grackle_field(HeI_density);
  rt_print_grackle_field(HeII_density);
  rt_print_grackle_field(HeIII_density);
  rt_print_grackle_field(e_density);
  rt_print_grackle_field(HM_density);
  rt_print_grackle_field(H2I_density);
  rt_print_grackle_field(H2II_density);
  rt_print_grackle_field(DI_density);
  rt_print_grackle_field(DII_density);
  rt_print_grackle_field(HDI_density);
  rt_print_grackle_field(metal_density);
  rt_print_grackle_field(x_velocity);
  rt_print_grackle_field(y_velocity);
  rt_print_grackle_field(z_velocity);
  rt_print_grackle_field(volumetric_heating_rate);
  rt_print_grackle_field(specific_heating_rate);
  rt_print_grackle_field(RT_HI_ionization_rate);
  rt_print_grackle_field(RT_HeI_ionization_rate);
  rt_print_grackle_field(RT_HeII_ionization_rate);
  rt_print_grackle_field(RT_H2_dissociation_rate);
  rt_print_grackle_field(RT_heating_rate);
}

#endif /* MY_GRACKLE_UTILS_H */
