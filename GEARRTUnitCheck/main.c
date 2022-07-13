/* ---------------------------------------------
 * In this example, we start with high internal
 * energies and a fully ionized gas, and just
 * let it cool without any RT.
 * --------------------------------------------- */

/* define these before including local headers like my_grackle_utils.h */
#define RT_NGROUPS 4
/* Grackle related macros */
#define FIELD_SIZE 1
#define GRIDDIM 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "parser.h"

/* #include <grackle.h> */
/*  */
#include "constants.h"
/* #include "cross_sections.h" */
/* #include "ionization_equilibrium.h" */
/* #include "mean_molecular_weight.h" */
/* #include "my_grackle_utils.h" */
/* #include "photon_interaction_rates.h" */
/*  */

/* Some global variables */
/* --------------------- */

/* Units to be used in the swift simulation */
double mass_units;
double time_units;
double length_units;
double density_units;
double velocity_units;
double temperature_units;
double energy_units;
double internal_energy_units;

double c_internal_units;

/* Radiation variables */
double *star_emission_rates;
double *photon_groups_Hz;

/* Other quantities. All assumed in SWIFT internal units. */
double particle_mass;
/* (Estimate of) minimal density in sim */
double density_min;
/* (Estimate of) maximal density in sim */
double density_max;
/* Average density in sim */
double density_average;
/* boxsize */
double boxsize;
/* (Estimate of) smoothing length */
double smoothing_length;
/* Number of gas particles in simulation */
long long npart;

/**
 * @brief Read in the parameters relevant for this check.
 *
 * @param params (return) swift_params struct to be filled
 * @param param_filenme filename to read in
 **/
void read_paramfile(struct swift_params *params, char *param_filename){
  message("Reading parameters from file '%s'", param_filename);
  parser_read_file(param_filename, params);
}

/**
 * @brief Read in the parameters relevant for this check.
 *
 * @param params swift_params struct to read from
 **/
void read_swift_params(struct swift_params *params){

  /* Read in data */
  mass_units = parser_get_param_double(params, "InternalUnitSystem:UnitMass_in_cgs");
  length_units = parser_get_param_double(params, "InternalUnitSystem:UnitLength_in_cgs");
  velocity_units = parser_get_param_double(params, "InternalUnitSystem:UnitVelocity_in_cgs");
  temperature_units = parser_get_param_double(params, "InternalUnitSystem:UnitTemp_in_cgs");

  double fc = parser_get_param_double(params, "GEARRT:f_reduce_c");
  c_internal_units = fc * const_speed_light_c / velocity_units;
  star_emission_rates = malloc(RT_NGROUPS * sizeof(double));
  double *photon_groups_read = malloc(RT_NGROUPS * sizeof(double));
  parser_get_param_double_array(params, "GEARRT:star_emission_rates_LSol", RT_NGROUPS, star_emission_rates);
  if (RT_NGROUPS == 0){
    error("Can't run RT with 0 photon groups. Modify RT_NGROUPS in this script.");
  } else if (RT_NGROUPS == 1){
    photon_groups_read[0] = 0.;
  } else{
    parser_get_param_double_array(params, "GEARRT:photon_groups_Hz", RT_NGROUPS-1, photon_groups_read);
  }
  photon_groups_Hz = malloc(RT_NGROUPS * sizeof(double));
  for (int i = 0; i < RT_NGROUPS-1; i++){
    photon_groups_Hz[i+1] = photon_groups_read[i];
  }
  photon_groups_Hz[0] = 0.;

  /* Convert units */
  time_units = length_units / velocity_units;
  density_units = mass_units / (length_units * length_units * length_units);
  internal_energy_units = velocity_units * velocity_units;
  energy_units = mass_units * internal_energy_units;

  /* Clean up */
  free(photon_groups_read);
}

/* Check that the value is a valid float.
 * Assume the argument given is positive and of type double. */
#define check_valid_float(v)                                            \
  ({                                                                    \
   if (v < 0.) {                                                        \
    fflush(stdout);                                                     \
    fprintf(stderr, "%s %s:%d " #v " invalid float; negative: %.6e\n",  \
        __FILE__, __FUNCTION__, __LINE__, v );                          \
    abort();                                                            \
   }                                                                    \
   if (v > FLT_MAX) {                                                   \
    fflush(stdout);                                                     \
    fprintf(stderr, "%s %s:%d " #v " invalid float; > FLT_MAX: %.6e\n", \
        __FILE__, __FUNCTION__, __LINE__, v);                           \
    abort();                                                            \
   }                                                                    \
   if ((v != 0.) && (v < FLT_MIN)) {                                    \
    fflush(stdout);                                                     \
    fprintf(stderr, "%s %s:%d " #v " invalid float; < FLT_MIN: %.6e\n", \
        __FILE__, __FUNCTION__, __LINE__, v);                           \
    abort();                                                            \
   }                                                                    \
   float vfloat = (float) v;                                            \
   if (isinf(vfloat) || isnan(vfloat)){                                 \
    fflush(stdout);                                                     \
    fprintf(stderr, "%s %s:%d " #v " invalid float; nan/inf %.6e\n",    \
        __FILE__, __FUNCTION__, __LINE__, v);                           \
    abort();                                                            \
   }                                                                    \
  })

/**
 * @brief Read in the parameters from the ICs relevant for this check.
 *
 * @param params swift_params struct to read from
 **/
void read_ic_params(struct swift_params *params){

  const double mass_units_ic = parser_get_param_double(params, "InternalUnitSystem:UnitMass_in_cgs");
  const double length_units_ic = parser_get_param_double(params, "InternalUnitSystem:UnitLength_in_cgs");
  /* const double velocity_units_ic = parser_get_param_double(params, "InternalUnitSystem:UnitVelocity_in_cgs"); */
  /* const double temperature_units_ic = parser_get_param_double(params, "InternalUnitSystem:UnitTemp_in_cgs"); */
  const double density_units_ic = mass_units_ic / (length_units_ic * length_units_ic * length_units_ic);

  const double particle_mass_ic = parser_get_param_double(params, "ParticleData:ParticleMass");
  const double av_density_ic = parser_get_param_double(params, "ParticleData:averageDensity");
  const double boxsize_ic = parser_get_param_double(params, "GlobalData:boxsize");
  npart = parser_get_param_longlong(params, "GlobalData:npart");

  /* Convert quantities from IC internal units to SWIFT internal units */
  message("Converting IC units to SWIFT run units");
  particle_mass = particle_mass_ic * mass_units_ic / mass_units;
  check_valid_float(particle_mass);

  density_average = av_density_ic * density_units_ic / density_units;
  check_valid_float(density_average);

  boxsize = boxsize_ic * length_units_ic / length_units;
  check_valid_float(boxsize);
  if (boxsize == 0.) error("Got 0. boxsize?");

  if (density_average == 0.){
    /* If density = 0 in parameter file, make an estimate yourself. */
    const double totmass = npart * particle_mass;
    density_average = totmass / (boxsize * boxsize * boxsize);
    check_valid_float(density_average);

    density_min = 1e-3 * density_average;
    check_valid_float(density_min);

    density_max = 1e3 * density_average;
    check_valid_float(density_max);
  } 
  else {
    const double min_density_ic = parser_get_param_double(params, "ParticleData:minDensity");
    density_min = min_density_ic * density_units_ic / density_units;
    check_valid_float(density_min);

    const double max_density_ic = parser_get_param_double(params, "ParticleData:maxDensity");
    density_max = max_density_ic * density_units_ic / density_units;
    check_valid_float(density_max);
  }

  const double sml_ic = parser_get_param_double(params, "ParticleData:smoothingLength");
  if (sml_ic == 0.){
    /* Try and estimate yourself */
    smoothing_length = 0.75 * boxsize / pow(npart, 0.3333333);
    check_valid_float(smoothing_length);
  }
  else {
    smoothing_length = sml_ic * length_units_ic / length_units;
    check_valid_float(smoothing_length);
  }
}


int main() {
  /* TODO: make this a cmdline arg? */
  /* char *swift_param_filename = "swift_parameters.yml"; */
  char *swift_param_filename = "ilievTest0part3.yml";
  char *sim_param_filename = "simulation_parameters.yml";

  /* Read the SWIFT parameter file, and the parameters */
  struct swift_params *swift_params =
      (struct swift_params *)malloc(sizeof(struct swift_params));
  read_paramfile(swift_params, swift_param_filename);
  read_swift_params(swift_params);

  /* Read the simulation data parameter file, and the parameters */
  struct swift_params *sim_params =
      (struct swift_params *)malloc(sizeof(struct swift_params));
  read_paramfile(sim_params, sim_param_filename);
  read_ic_params(sim_params);

  /* Clean up after yourself */
  free(swift_params);
  free(sim_params);
  free(star_emission_rates);
  free(photon_groups_Hz);
  return 0;
}
