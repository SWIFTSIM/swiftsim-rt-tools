#include "constants.h"
#include "error.h"
#include "units.h"
#include "simulation_params.h"
#include "validity_check_macros.h"

extern int warnings;


/**
 * @brief Read in the parameters relevant for this check from
 * the parameter file into the swift_parames struct.
 *
 * @param params (return) swift_params struct to be filled
 * @param param_filenme filename to read in
 **/
void params_read_paramfile(struct swift_params *params, char *param_filename) {
  message("Reading parameters from file '%s'", param_filename);
  parser_read_file(param_filename, params);
}

/**
 * @brief Read in the parameters relevant for this check from the
 * swift_params struct.
 *
 * @param params swift_params struct to read from
 * @param units units to be used in swift internally during the run.
 * @param simulation_params struct containing parameters we'll actually need in this check
 **/
void params_read_swift_params(struct swift_params *params, struct units* units, struct simulation_params* simulation_params) {

  /* Read in data */
  double mass_units =
      parser_get_param_double(params, "InternalUnitSystem:UnitMass_in_cgs");
  double length_units =
      parser_get_param_double(params, "InternalUnitSystem:UnitLength_in_cgs");
  double velocity_units =
      parser_get_param_double(params, "InternalUnitSystem:UnitVelocity_in_cgs");
  double temperature_units =
      parser_get_param_double(params, "InternalUnitSystem:UnitTemp_in_cgs");

  float fc = parser_get_param_float(params, "GEARRT:f_reduce_c");
  double c_reduced = fc * const_speed_light_c / velocity_units;
  double* photon_groups_Hz = malloc(RT_NGROUPS * sizeof(double));
  if (RT_NGROUPS == 0) {
    error(
        "Can't run RT with 0 photon groups. Modify RT_NGROUPS in this script.");
  } else {
    parser_get_param_double_array(params, "GEARRT:photon_groups_Hz", RT_NGROUPS,
                                  photon_groups_Hz);
  }

  int use_const_emission_rates = 0;

  char stellar_model_str[80];
  parser_get_param_string(params, "GEARRT:stellar_luminosity_model",
                          stellar_model_str);
  if (strcmp(stellar_model_str, "const") == 0)
    use_const_emission_rates = 1;

  double *star_emission_rates = malloc(RT_NGROUPS * sizeof(double));

  if (use_const_emission_rates) {
    parser_get_param_double_array(params,
                                  "GEARRT:const_stellar_luminosities_LSol",
                                  RT_NGROUPS, star_emission_rates);

  } else {
    error("This check isn't set up to run without constant stellar emission "
          "rates (yet)");
  }


  /* Now copy the data into the place they belong. */

  units->mass_units = mass_units;
  units->length_units = length_units;
  units->velocity_units = velocity_units;
  units->temperature_units = temperature_units;

  simulation_params->fc = fc;
  simulation_params->c_reduced = c_reduced;
  for (int g = 0; g < RT_NGROUPS; g++){
    simulation_params->photon_groups_Hz[g] = photon_groups_Hz[g];
    simulation_params->star_emission_rates[g] = star_emission_rates[g];
  }
  simulation_params->use_const_emission_rates = use_const_emission_rates;


}


/**
 * @brief Read in the parameters from the ICs relevant for this check.
 *
 * @param params swift_params struct to read from
 * @param units units to be used in swift internally during the run.
 * @param simulation_params struct containing parameters we'll actually need in this check
 * @param verbose How verbose to be.
 **/
void params_read_ic_params(struct swift_params *params, struct units* units, struct simulation_params* simulation_params, int verbose) {

  const double mass_units_ic =
      parser_get_param_double(params, "InternalUnitSystem:UnitMass_in_cgs");
  const double length_units_ic =
      parser_get_param_double(params, "InternalUnitSystem:UnitLength_in_cgs");
  const double velocity_units_ic =
      parser_get_param_double(params, "InternalUnitSystem:UnitVelocity_in_cgs");
  /* const double temperature_units_ic = parser_get_param_double(params,
   * "InternalUnitSystem:UnitTemp_in_cgs"); */
  const double density_units_ic =
      mass_units_ic / (length_units_ic * length_units_ic * length_units_ic);
  const double energy_units_ic =
      mass_units_ic * velocity_units_ic * velocity_units_ic;

  const float particle_mass_ic =
      parser_get_param_float(params, "ParticleData:ParticleMass");
  const float av_density_ic =
      parser_get_param_float(params, "ParticleData:averageDensity");
  const float av_radiation_ic =
      parser_get_param_float(params, "ParticleData:averageRadiationEnergy");
  const float boxsize_ic = parser_get_param_float(params, "GlobalData:boxsize");
  const long long npart = parser_get_param_longlong(params, "GlobalData:npart");

  /* Convert quantities from IC internal units to SWIFT internal units */
  if (verbose) message("Converting IC units to SWIFT run units");

  const float particle_mass = particle_mass_ic * mass_units_ic / units->mass_units;
  check_valid_float(particle_mass, 0);

  float density_average = av_density_ic * density_units_ic / units->density_units;
  check_valid_float(density_average, 0);

  if (density_average > 0.)
    check_valid_float(density_average, 1);

  const float rad_energy_av = av_radiation_ic * energy_units_ic / units->energy_units;
  if (rad_energy_av != 0.)
    check_valid_float(rad_energy_av, 0);
  if (rad_energy_av < 0.)
    error("Average radiation energy = %.3g < 0", rad_energy_av);

  const float boxsize = boxsize_ic * length_units_ic / units->length_units;
  check_valid_float(boxsize, 1);
  if (boxsize == 0.)
    error("Got 0. boxsize?");

  float density_min = 0.f;
  float density_max = 0.f;

  if (density_average == 0.) {
    /* If density = 0 in parameter file, make an estimate yourself. */
    if (verbose) message("Got zero average density from IC paramfile, estimating myself.");

    const double totmass = npart * particle_mass;
    density_average = totmass / (boxsize * boxsize * boxsize);
    check_valid_float(density_average, 1);

    density_min = 1e-3 * density_average;
    check_valid_float(density_min, 1);

    density_max = 1e3 * density_average;
    check_valid_float(density_max, 1);

  } else {

    check_valid_float(density_average, 1);

    const float min_density_ic =
        parser_get_param_float(params, "ParticleData:minDensity");
    if (min_density_ic > av_density_ic)
      error("density min > average? %.6e %.6e", min_density_ic, av_density_ic);

    const float max_density_ic =
        parser_get_param_float(params, "ParticleData:maxDensity");
    if (max_density_ic < av_density_ic)
      error("density max < average? %.6e %.6e", min_density_ic, av_density_ic);

    density_min = min_density_ic * density_units_ic / units->density_units;
    check_valid_float(density_min, 1);

    density_max = max_density_ic * density_units_ic / units->density_units;
    check_valid_float(density_max, 1);
  }

  float smoothing_length = 0.f;
  const float sml_ic =
      parser_get_param_float(params, "ParticleData:smoothingLength");
  if (sml_ic == 0.) {
    /* Try and estimate yourself */
    if (verbose) message("Got zero smoothing length from IC paramfile, estimating myself.");
    smoothing_length = 0.75 * boxsize / pow(npart, 0.3333333);
    check_valid_float(smoothing_length, 0);
  } else {
    smoothing_length = sml_ic * length_units_ic / units->length_units;
    check_valid_float(smoothing_length, 0);
  }

  float rad_energy_min = 0.f;
  float rad_energy_max = 0.f;
  if (rad_energy_av > 0.) {
    const double rad_energy_max_ic =
        parser_get_param_float(params, "ParticleData:maxRadiationEnergy");
    rad_energy_max = rad_energy_max_ic * energy_units_ic / units->energy_units;
    const double rad_energy_min_ic =
        parser_get_param_float(params, "ParticleData:minRadiationEnergy");
    rad_energy_min = rad_energy_min_ic * energy_units_ic / units->energy_units;
  }


  /* ------------------------------------------------ */
  /* Finally, store the parameters where they belong. */
  /* ------------------------------------------------ */

  simulation_params->particle_mass = particle_mass;
  simulation_params->density_average = density_average;
  simulation_params->density_min = density_min;
  simulation_params->density_max = density_max;
  simulation_params->boxsize = boxsize;
  simulation_params->smoothing_length = boxsize;

  simulation_params->rad_energy_min = rad_energy_min;
  simulation_params->rad_energy_max = rad_energy_max;
  simulation_params->rad_energy_av = rad_energy_av;

  simulation_params->npart = npart;


}


