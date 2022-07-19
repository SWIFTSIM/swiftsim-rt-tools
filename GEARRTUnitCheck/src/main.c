/* -------------------------------------------------
 * A suite of tests both with and without grackle
 * to check whether your choice of units may produce
 * bad results.
 * ------------------------------------------------- */

/* FPE stuff. Link with -lm */
/* #define _GNU_SOURCE  */
/* #include <fenv.h>  */

/* define these before including local headers like my_grackle_utils.h */
#define RT_NGROUPS 4
/* Grackle related macros */
#define FIELD_SIZE 1
#include "my_grackle_utils.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "grackle_cooling_test.h"
#include "grackle_heating_test.h"
#include "parser.h"

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

/* Radiation variables */
float c_reduced; /* in internal units */
double *star_emission_rates;
double *photon_groups_Hz;

/* Other quantities. All assumed in SWIFT internal units. */
float particle_mass;
/* (Estimate of) minimal density in sim */
float density_min;
/* (Estimate of) maximal density in sim */
float density_max;
/* Average density in sim */
float density_average;
/* boxsize */
float boxsize;
/* (Estimate of) smoothing length */
float smoothing_length;
/* Number of gas particles in simulation */
long long npart;

int warnings = 0;

/*! Check that the value is a valid float.
 * Assume the argument given is positive and of type double. */
#define check_valid_float(v)                                                   \
  ({                                                                           \
    if (v < 0.) {                                                              \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid float; negative: %.6e\n",      \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if (v > FLT_MAX) {                                                         \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid float; > FLT_MAX: %.6e\n",     \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if ((v != 0.) && (v < FLT_MIN)) {                                          \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid float; < FLT_MIN: %.6e\n",     \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    float vfloat = (float)v;                                                   \
    if (isinf(vfloat) || isnan(vfloat)) {                                      \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid float; nan/inf %.6e\n",        \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if (v != 0. && (fabs(v) > 1e30 || fabs(v) < 1e-30)) {                      \
      fflush(stdout);                                                          \
      fprintf(stdout, "WARNING: %s:%s:%d: " #v " has large exponent: %.6e\n",  \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      warnings++;                                                              \
    }                                                                          \
  })

#define FLOAT_TOLERANCE 1e-5
#define check_floats_equal(a, b)                                               \
  ({                                                                           \
    if (a == 0. && b == 0.) {                                                  \
    } else if (a == 0. && fabsf(b) > FLOAT_TOLERANCE) {                        \
      error(#a " and " #b " are not equal: %.6e %.6e", a, b);                  \
    } else if (b == 0. && fabsf(a) > FLOAT_TOLERANCE) {                        \
      error(#a " and " #b " are not equal: %.6e %.6e", a, b);                  \
    } else {                                                                   \
      if (1.f - fabsf(a / b) > FLOAT_TOLERANCE) {                              \
        error(#a " and " #b " are not equal: %.6e %.6e a/b=%.3e", a, b,        \
              a / b);                                                          \
      }                                                                        \
    }                                                                          \
  })

/*! Check that the value is a valid double.
 * Assume the argument given is positive and of type double. */
#define check_valid_double(v)                                                  \
  ({                                                                           \
    if (v < 0.) {                                                              \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid double; negative: %.6e\n",     \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if (isinf(v) || isnan(v)) {                                                \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid double; nan/inf %.6e\n",       \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if (v > DBL_MAX) {                                                         \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid double; > DBL_MAX: %.6e\n",    \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if ((v != 0.) && (v < DBL_MIN)) {                                          \
      fflush(stdout);                                                          \
      fprintf(stderr, "%s:%s:%d: " #v " invalid float; < DBL_MIN: %.6e\n",     \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      abort();                                                                 \
    }                                                                          \
    if (v != 0. && (fabs(v) > 1e290 || fabs(v) < 1e-290)) {                    \
      fflush(stdout);                                                          \
      fprintf(stdout, "WARNING: %s:%s:%d: " #v " has large exponent: %.6e\n",  \
              __FILE__, __FUNCTION__, __LINE__, v);                            \
      warnings++;                                                              \
    }                                                                          \
  })

#define DOUBLE_TOLERANCE 1e-15
#define check_doubles_equal(a, b)                                               \
  ({                                                                           \
    if (a == 0. && b == 0.) {                                                  \
    } else if (a == 0. && fabs(b) > DOUBLE_TOLERANCE) {                        \
      error(#a " and " #b " are not equal: %.6e %.6e", a, b);                  \
    } else if (b == 0. && fabs(a) > DOUBLE_TOLERANCE) {                        \
      error(#a " and " #b " are not equal: %.6e %.6e", a, b);                  \
    } else {                                                                   \
      if (1.f - fabs(a / b) > DOUBLE_TOLERANCE) {                              \
        error(#a " and " #b " are not equal: %.6e %.6e a/b=%.3e", a, b,        \
              a / b);                                                          \
      }                                                                        \
    }                                                                          \
  })

/**
 * @brief Read in the parameters relevant for this check.
 *
 * @param params (return) swift_params struct to be filled
 * @param param_filenme filename to read in
 **/
void read_paramfile(struct swift_params *params, char *param_filename) {
  message("Reading parameters from file '%s'", param_filename);
  parser_read_file(param_filename, params);
}

/**
 * @brief Read in the parameters relevant for this check.
 *
 * @param params swift_params struct to read from
 **/
void read_swift_params(struct swift_params *params) {

  /* Read in data */
  mass_units =
      parser_get_param_double(params, "InternalUnitSystem:UnitMass_in_cgs");
  length_units =
      parser_get_param_double(params, "InternalUnitSystem:UnitLength_in_cgs");
  velocity_units =
      parser_get_param_double(params, "InternalUnitSystem:UnitVelocity_in_cgs");
  temperature_units =
      parser_get_param_double(params, "InternalUnitSystem:UnitTemp_in_cgs");

  float fc = parser_get_param_float(params, "GEARRT:f_reduce_c");
  c_reduced = fc * const_speed_light_c / velocity_units;
  star_emission_rates = malloc(RT_NGROUPS * sizeof(double));
  double *photon_groups_read = malloc(RT_NGROUPS * sizeof(double));
  parser_get_param_double_array(params, "GEARRT:star_emission_rates_LSol",
                                RT_NGROUPS, star_emission_rates);
  if (RT_NGROUPS == 0) {
    error(
        "Can't run RT with 0 photon groups. Modify RT_NGROUPS in this script.");
  } else if (RT_NGROUPS == 1) {
    photon_groups_read[0] = 0.;
  } else {
    parser_get_param_double_array(params, "GEARRT:photon_groups_Hz",
                                  RT_NGROUPS - 1, photon_groups_read);
  }
  photon_groups_Hz = malloc(RT_NGROUPS * sizeof(double));
  for (int i = 0; i < RT_NGROUPS - 1; i++) {
    photon_groups_Hz[i + 1] = photon_groups_read[i];
  }
  photon_groups_Hz[0] = 0.;

  /* Convert units */
  time_units = length_units / velocity_units;
  check_valid_float(time_units);
  density_units = mass_units / (length_units * length_units * length_units);
  check_valid_float(density_units);
  internal_energy_units = velocity_units * velocity_units;
  check_valid_float(internal_energy_units);
  energy_units = mass_units * internal_energy_units;
  check_valid_float(energy_units);

  /* Clean up */
  free(photon_groups_read);
}
/**
 * @brief Read in the parameters from the ICs relevant for this check.
 *
 * @param params swift_params struct to read from
 **/
void read_ic_params(struct swift_params *params) {

  const double mass_units_ic =
      parser_get_param_double(params, "InternalUnitSystem:UnitMass_in_cgs");
  const double length_units_ic =
      parser_get_param_double(params, "InternalUnitSystem:UnitLength_in_cgs");
  /* const double velocity_units_ic = parser_get_param_double(params,
   * "InternalUnitSystem:UnitVelocity_in_cgs"); */
  /* const double temperature_units_ic = parser_get_param_double(params,
   * "InternalUnitSystem:UnitTemp_in_cgs"); */
  const double density_units_ic =
      mass_units_ic / (length_units_ic * length_units_ic * length_units_ic);

  const float particle_mass_ic =
      parser_get_param_float(params, "ParticleData:ParticleMass");
  const float av_density_ic =
      parser_get_param_float(params, "ParticleData:averageDensity");
  const float boxsize_ic = parser_get_param_float(params, "GlobalData:boxsize");
  npart = parser_get_param_longlong(params, "GlobalData:npart");

  /* Convert quantities from IC internal units to SWIFT internal units */
  message("Converting IC units to SWIFT run units");
  particle_mass = particle_mass_ic * mass_units_ic / mass_units;
  check_valid_float(particle_mass);

  density_average = av_density_ic * density_units_ic / density_units;
  check_valid_float(density_average);

  boxsize = boxsize_ic * length_units_ic / length_units;
  check_valid_float(boxsize);
  if (boxsize == 0.)
    error("Got 0. boxsize?");

  if (density_average == 0.) {
    /* If density = 0 in parameter file, make an estimate yourself. */
    const float totmass = npart * particle_mass;
    density_average = totmass / (boxsize * boxsize * boxsize);
    check_valid_float(density_average);

    density_min = 1e-3 * density_average;
    check_valid_float(density_min);

    density_max = 1e3 * density_average;
    check_valid_float(density_max);
  } else {
    const float min_density_ic =
        parser_get_param_float(params, "ParticleData:minDensity");
    if (min_density_ic > av_density_ic)
      error("density min > average? %.6e %.6e", min_density_ic, av_density_ic);

    const float max_density_ic =
        parser_get_param_float(params, "ParticleData:maxDensity");
    if (max_density_ic < av_density_ic)
      error("density max < average? %.6e %.6e", min_density_ic, av_density_ic);

    density_min = min_density_ic * density_units_ic / density_units;
    check_valid_float(density_min);

    density_max = max_density_ic * density_units_ic / density_units;
    check_valid_float(density_max);
  }

  const float sml_ic =
      parser_get_param_float(params, "ParticleData:smoothingLength");
  if (sml_ic == 0.) {
    /* Try and estimate yourself */
    smoothing_length = 0.75 * boxsize / pow(npart, 0.3333333);
    check_valid_float(smoothing_length);
  } else {
    smoothing_length = sml_ic * length_units_ic / length_units;
    check_valid_float(smoothing_length);
  }
}

/**
 * @brief print out the used parameters for a visual inspection
 **/

void print_params() {

  message("Units: [in cgs]");
  message("%22s: %.6e", "mass units", mass_units);
  message("%22s: %.6e", "length units", length_units);
  message("%22s: %.6e", "time units", time_units);
  message("%22s: %.6e", "density units", density_units);
  message("%22s: %.6e", "velocity units", velocity_units);
  message("%22s: %.6e", "temperature units", temperature_units);
  message("%22s: %.6e", "energy units", energy_units);
  message("%22s: %.6e", "internal energy units", internal_energy_units);

  message("");
  message("Variables [internal units]");
  message("%22s: %.6e", "particle mass", particle_mass);
  message("%22s: %.6e", "density_average", density_average);
  message("%22s: %.6e", "density_min", density_min);
  message("%22s: %.6e", "density_max", density_max);
  message("%22s: %.6e", "boxsize", boxsize);
  message("%22s: %.6e", "smoothing length", smoothing_length);
  message("%22s: %.6e", "approx dt [internal units]",
          smoothing_length / c_reduced);
  message("%22s: %.6e", "approx dt [s]             ",
          smoothing_length / c_reduced * time_units);
  message("%22s: %.6e", "approx dt [kyr]           ",
          smoothing_length / c_reduced * time_units / const_yr * 1e-3);
}

/**
 * @brief Check whether other gas quantities derived from the
 * initial conditions are within an acceptable range
 *
 * @param density gas density to use
 * @param name name of the test case. NO SPACES.
 * @param float T temperature to deal with
 * @param verbose are we talkative?
 **/
void check_gas_quantities(float density, char *name, float T, int verbose) {

  /* assume mean molecular weight of 1 for this test. While that isn't correct,
   * it should do the trick for the purpose of this test. */
  message("Checking gas quantities for T=%.1f case=%s", T, name);

  const float gamma = const_adiabatic_index;
  const float gamma_minus_one = gamma - 1.f;

  /* Get and check internal energy */
  const float mu = 1.;
  const float internal_energy_cgs =
      const_kboltz * T / (gamma_minus_one * mu * const_mh);
  const float internal_energy = internal_energy_cgs / internal_energy_units;
  check_valid_float((double)internal_energy);

  /* Get and check other quantities from internal energy */
  const float pressure = gamma_minus_one * internal_energy * density;
  check_valid_float((double)pressure);

  const float entropy = pressure * pow(density, -gamma);
  check_valid_float((double)entropy);

  const float soundspeed = sqrtf(internal_energy * gamma * gamma_minus_one);
  check_valid_float((double)soundspeed);

  /* Get and check quantities using entropy now */
  const float internal_energy_from_entropy =
      entropy * powf(density, gamma_minus_one) / gamma_minus_one;
  check_valid_float((double)internal_energy_from_entropy);

  const float pressure_from_entropy = entropy * powf(density, gamma);
  check_valid_float((double)pressure_from_entropy);

  const float entropy_from_internal_energy =
      gamma_minus_one * internal_energy * powf(density, -gamma_minus_one);
  check_valid_float(entropy_from_internal_energy);

  const float soundspeed_from_entropy =
      sqrtf(gamma * powf(density, gamma_minus_one) * entropy);
  check_valid_float(soundspeed_from_entropy);

  const float internal_energy_from_pressure =
      (pressure_from_entropy / density) / gamma_minus_one;
  check_valid_float((double)internal_energy_from_pressure);

  const float soundspeed_from_pressure = sqrtf(gamma * pressure / density);
  check_valid_float((double)soundspeed_from_pressure);

  /* Compare obtained values */
  check_floats_equal(pressure, pressure_from_entropy);
  check_floats_equal(entropy, entropy_from_internal_energy);
  check_floats_equal(internal_energy, internal_energy_from_entropy);
  check_floats_equal(internal_energy, internal_energy_from_pressure);
  check_floats_equal(soundspeed, soundspeed_from_entropy);
  check_floats_equal(soundspeed, soundspeed_from_pressure);

  /* Try and obtain converved fluid quantities using the local soundspeed
   * as an estimate for fluid velocity */

  const double momentum = particle_mass * soundspeed;
  check_valid_float(momentum);

  const double total_energy =
      particle_mass * (internal_energy + 0.5 * soundspeed * soundspeed);
  check_valid_float(total_energy);

  if (verbose) {
    message("rho=%.3e u=%.3e P=%.3e A=%.3e cs=%.3e | mass=%.3e momentum=%.3e "
            "E=%.3e",
            density, internal_energy, pressure, entropy, soundspeed,
            particle_mass, momentum, total_energy);
  }

  /* Check that the number densities are well behaved.
   * The assumption here is that we never use the number densities
   * as floats. Assume pure hydrogen gas. */

  const double n = density / (const_mh / mass_units);
  check_valid_double(n);
  const double n_cgs = density * density_units / const_mh;
  check_valid_double(n_cgs);
  /* Make sure the number density is within the limit of what grackle can do */
  if (n_cgs < 1e-10)
    error("case=%s density=%.3e gives number density=%.3e [cm^-3] which is "
          "below lower limit of 1e-10",
          name, density, n_cgs);
  if (n_cgs > 1e16)
    error("case=%s density=%.3e gives number density=%.3e [cm^-3] which is "
          "above upper limit of 1e16",
          name, density, n_cgs);

  /* Check for heating rates divided by number density, which is a parameter 
   * that grackle requires. Pretend we have hydrogen only gas, and only 1 
   * radiation group.  Then try with "fully neutral" and "fully ionized" gas,
   * where we assume that "fully ionized" is taken to mean that the number
   * density of neutral hydrogen is multiplied by a factor of 1e-9. 
   * The cross sections assume a blackbody spectrum of T = 1e5K.
   * */

  const double cse = 1.096971e-18;
  const double csn = 1.630511e-18;
  const double mean_photon_energy = 4.744e-11;
  const double fixed_luminosity_cgs = 4.774e+01; /* erg/s/cm^2 */
  const double E_ion = 2.179e-11; /* erg, ionization energy of hydrogen */

  const double fixed_radiation_density_cgs = fixed_luminosity_cgs / const_speed_light_c;
  const double Eic = fixed_radiation_density_cgs * const_speed_light_c;
  const double Nic = Eic / mean_photon_energy;

  double ionization_factor[2] = {1., 1.e-9};
  double heating_rates[2] = {0., 0.};

  for (int i = 0; i < 2; i++) {
    const double f = ionization_factor[i];
    const double HI_heating_rate = (cse * mean_photon_energy - E_ion * csn) * Nic * (n_cgs * f);
    const double HI_ionization_rate = csn * (n_cgs * f) * Nic;

    check_valid_double(HI_heating_rate);
    check_valid_double(HI_ionization_rate);
    /* heating rate is expected to be in erg/s/cm^3 / nHI*/
    const double HI_heating_rate_for_grackle = HI_heating_rate / (n * f);
    /* ion. rate is expected to be in 1/time_units, so divide by (1/time_units) */
    const double HI_ionization_rate_for_grackle = HI_ionization_rate * time_units;

    if (verbose) {
      message("Heating rate: %.4e Ionization rate: %.4e", 
          HI_heating_rate, HI_ionization_rate);
      message("Heating rate GR: %.4e Ionization rate GR: %.4e", 
          HI_heating_rate_for_grackle, HI_ionization_rate_for_grackle);
    }

    check_valid_double(HI_heating_rate_for_grackle);
    check_valid_double(HI_ionization_rate_for_grackle);
    heating_rates[i] = HI_heating_rate_for_grackle;
  }
  check_doubles_equal(heating_rates[0], heating_rates[1]);
}

int main(void) {

  /* FPE's, here we come! */
  /* Note: Grackle may have severe issues with these... And that's apparently
   * ok? I gave up. */
  /* feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW); */

  /* ----------------------------------------- */
  /* First things first: Read in required data */
  /* ----------------------------------------- */

  /* TODO: make this a cmdline arg? */
  /* char *swift_param_filename = "swift_parameters.yml"; */
  char *swift_param_filename = "ilievTest0part3.yml";
  char *sim_param_filename = "simulation_parameters.yml";

  /* Print a lot of information to the screen? */
  int verbose = 0;

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

  /* Print out the units and parameters for a visual inspection */
  print_params();

  /* ------------------------*/
  /* Prepare to run examples */
  /* ------------------------*/
  /* If the min/max densities are too close to the average, re-size them */

  if (fabs(1. - density_min / density_average) < 0.05) {
    message("density_min too close to average. Resizing %.3e -> %.3e",
            density_min, 0.8 * density_average);
    density_min = 0.8 * density_average;
    check_valid_float(density_min);
  }
  if (fabs(density_max / density_average - 1.) < 0.05) {
    message("density_max too close to average. Resizing %.3e -> %.3e",
            density_max, 1.2 * density_average);
    density_max = 1.2 * density_average;
    check_valid_float(density_max);
  }

  float dens_arr[3] = {density_average, density_min, density_max};
  char *dens_names[3] = {"density_average", "density_min", "density_max"};

  /* Run Gas Quantities Checks */
  /* ------------------------- */
  for (int d = 0; d < 3; d++) {
    float rho = dens_arr[d];
    char *name = dens_names[d];
    check_gas_quantities(rho, name, /*T=*/10., verbose);
    check_gas_quantities(rho, name, /*T=*/100., verbose);
    check_gas_quantities(rho, name, /*T=*/1000., verbose);
    check_gas_quantities(rho, name, /*T=*/10000., verbose);
    check_gas_quantities(rho, name, /*T=*/100000., verbose);
    check_gas_quantities(rho, name, /*T=*/1000000., verbose);
    check_gas_quantities(rho, name, /*T=*/10000000., verbose);
  }

  /* Run Grackle cooling test */
  /* ------------------------ */
  for (int d = 0; d < 3; d++) {
    float rho = dens_arr[d];
    char *name = dens_names[d];
    run_grackle_cooling_test(rho, name, mass_units, length_units, time_units,
                             density_units, velocity_units,
                             internal_energy_units, verbose);
  }

  /* Run Grackle heating test */
  /* ------------------------ */
  for (int d = 0; d < 3; d++) {
    float rho = dens_arr[d];
    char *name = dens_names[d];
    run_grackle_heating_test(rho, name, mass_units, length_units, time_units,
                             density_units, velocity_units,
                             internal_energy_units, verbose);
  }

  /* TODO: Luminosities check */

  /* Clean up after yourself */
  free(swift_params);
  free(sim_params);
  free(star_emission_rates);
  free(photon_groups_Hz);

  message("Check completed with %d warning(s).", warnings);
  message("Bye, and good luck.");
  return 0;
}
