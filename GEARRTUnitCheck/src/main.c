/* -------------------------------------------------
 * A suite of tests both with and without grackle
 * to check whether your choice of units may produce
 * bad results.
 * ------------------------------------------------- */

/* FPE stuff. Link with -lm */
/* #define _GNU_SOURCE  */
/* #include <fenv.h>  */

/* define RT_NGROUPS before including local headers like my_grackle_utils.h */
#include "rt_ngroups.h"

/* Grackle related macros */
#define FIELD_SIZE 1
/* don't raise warnings about deprecated grackle functions */
#define OMIT_LEGACY_INTERNAL_GRACKLE_FUNC
#include "my_grackle_utils.h"

#include <float.h>
#include <math.h>
/* #include <stdio.h> */
/* #include <stdlib.h> */
#include <string.h>

#include "constants.h"
#include "conversions.h"
#include "grackle_cooling_test.h"
#include "grackle_heating_test.h"
#include "ionization_equilibrium.h"
#include "mean_molecular_weight.h"
#include "validity_check_macros.h"
#include "units.h"
#include "read_params.h"

/* Some global variables */
/* --------------------- */

/* Units to be used in the swift simulation */

/* [> Radiation variables <] */
/* float c_reduced = 0.; [> in internal units <] */
/* double *star_emission_rates = NULL; */
/* double *photon_groups_Hz = NULL; */
/* int use_const_emission_rates = 0.; */
/*  */
/* [> Other quantities. All assumed in SWIFT internal units. <] */
/* float particle_mass = 0.f; */
/* [> (Estimate of) minimal density in sim <] */
/* float density_min = 0.f; */
/* [> (Estimate of) maximal density in sim <] */
/* float density_max = 0.f; */
/* [> Average density in sim <] */
/* float density_average = 0.f; */
/* [> boxsize <] */
/* float boxsize = 0.f; */
/* [> (Estimate of) smoothing length <] */
/* float smoothing_length = 0.f; */
/* [> Minimal radiation energy (not density) in ICs <] */
/* float rad_energy_min = 0.f; */
/* [> Maximal radiation energy (not density) in ICs <] */
/* float rad_energy_max = 0.f; */
/* [> Average radiation energy (not density) in ICs <] */
/* float rad_energy_av = 0.f; */
/* [> Number of gas particles in simulation <] */
/* long long npart = 0.l; */
/*  */
int warnings = 0;

/**
 * @brief Check whether other gas quantities derived from the
 * initial conditions are within an acceptable range
 *
 * @param density gas density to use
 * @param name name of the test case. NO SPACES.
 * @param T temperature to deal with
 * @param params simulation parameters.
 * @param units internal units used for the simulation
 * @param verbose are we talkative?
 **/
void check_gas_quantities(float density, char *name, float T, const struct simulation_params* params, const struct units* units, int verbose) {

  /* assume mean molecular weight of 1 for this test. While that isn't correct,
   * it should do the trick for the purpose of this test. */
  message("Checking T=%.1e case=%s", T, name);

  const float gamma = const_adiabatic_index;
  const float gamma_minus_one = gamma - 1.f;

  /* Get and check internal energy */
  const float mu = 1.;
  const float internal_energy_cgs =
      const_kboltz * T / (gamma_minus_one * mu * const_mh);
  const float internal_energy = internal_energy_cgs / units->internal_energy_units;
  check_valid_float((double)internal_energy, 1);

  /* Get and check other quantities from internal energy */
  const float pressure = gamma_minus_one * internal_energy * density;
  check_valid_float((double)pressure, 1);

  const float entropy = pressure * pow(density, -gamma);
  check_valid_float((double)entropy, 1);

  const float soundspeed = sqrtf(internal_energy * gamma * gamma_minus_one);
  check_valid_float((double)soundspeed, 1);

  /* Get and check quantities using entropy now */
  const float internal_energy_from_entropy =
      entropy * powf(density, gamma_minus_one) / gamma_minus_one;
  check_valid_float((double)internal_energy_from_entropy, 1);

  const float pressure_from_entropy = entropy * powf(density, gamma);
  check_valid_float((double)pressure_from_entropy, 1);

  const float entropy_from_internal_energy =
      gamma_minus_one * internal_energy * powf(density, -gamma_minus_one);
  check_valid_float(entropy_from_internal_energy, 1);

  const float soundspeed_from_entropy =
      sqrtf(gamma * powf(density, gamma_minus_one) * entropy);
  check_valid_float(soundspeed_from_entropy, 1);

  const float internal_energy_from_pressure =
      (pressure_from_entropy / density) / gamma_minus_one;
  check_valid_float((double)internal_energy_from_pressure, 1);

  const float soundspeed_from_pressure = sqrtf(gamma * pressure / density);
  check_valid_float((double)soundspeed_from_pressure, 1);

  /* Compare obtained values */
  check_floats_equal(pressure, pressure_from_entropy);
  check_floats_equal(entropy, entropy_from_internal_energy);
  check_floats_equal(internal_energy, internal_energy_from_entropy);
  check_floats_equal(internal_energy, internal_energy_from_pressure);
  check_floats_equal(soundspeed, soundspeed_from_entropy);
  check_floats_equal(soundspeed, soundspeed_from_pressure);

  /* Try and obtain converved fluid quantities using the local soundspeed
   * as an estimate for fluid velocity */

  const double momentum = params->particle_mass * soundspeed;
  check_valid_float(momentum, 1);

  const double total_energy =
      params->particle_mass * (internal_energy + 0.5 * soundspeed * soundspeed);
  check_valid_float(total_energy, 1);

  if (verbose) {
    message("rho=%.3e u=%.3e P=%.3e A=%.3e cs=%.3e | mass=%.3e momentum=%.3e "
            "E=%.3e",
            density, internal_energy, pressure, entropy, soundspeed,
            params->particle_mass, momentum, total_energy);
  }
}

/**
 * @brief attempt to check whether quantities used by
 * grackle internally have valid values.
 *
 * @param density gas density to use
 * @param radiation_energy_density radiation energy density to use
 * @param name name of the test case. NO SPACES.
 * @param T temperature to deal with
 * @param units the internal units to be used in the simulation
 * @param verbose are we talkative?
 **/
void check_grackle_internals(float density, float radiation_energy_density,
                             char *name, float T, const struct units* units, int verbose) {

  message("checking %s, T=%.1e", name, T);

  /* Check that the total number density is within the valid limits for
   * grackle to handle. */
  const double n_cgs = density * units->density_units / const_mh;
  check_valid_double(n_cgs, 0);
  /* Make sure the number density is within the limit of what grackle can do */
  if (n_cgs < 1e-10)
    error("case=%s density=%.3e gives number density=%.3e [cm^-3] which is "
          "below lower limit of 1e-10",
          name, density, n_cgs);
  if (n_cgs > 1e16)
    error("case=%s density=%.3e gives number density=%.3e [cm^-3] which is "
          "above upper limit of 1e16",
          name, density, n_cgs);

  /* Prepare some grackle internals for later. See cool1d_multi_g.F
   * in grackle's source files for reference. */
  const double dom = units->density_units / const_mh;
  check_valid_double(dom, 0);
  const double tbase1 = units->time_units;
  check_valid_double(tbase1, 0);
  const double xbase1 = units->length_units;
  check_valid_double(xbase1, 0);
  const double dbase1 = units->density_units;
  check_valid_double(dbase1, 0);
  const double coolunit = xbase1 * xbase1 * const_mh * const_mh /
                          (tbase1 * tbase1 * tbase1 * dbase1);
  check_valid_double(coolunit, 0);

  /* Test for different gas compositions. */
  double hydrogen_mass_fractions[3] = {1., 0.75, 0.25};
  /* Store heating rates for each composition. Grackle wants them in
   * units of erg/s/cm^2 / nHI_cgs, so for each mass fraction of
   * hydrogen, the heating rates that we provide to grackle should be
   * identical. Check whether that is the case after the computation
   * is done. */
  double heating_rates[3] = {0., 0., 0.};

  for (int x = 0; x < 3; x++) {
    /* Get and check that the number densities are well behaved.
     * The assumption here is that we never use the number densities
     * as floats. Assume pure hydrogen gas. */

    double X = hydrogen_mass_fractions[x];
    double nH = X * density / (const_mh / units->mass_units);
    check_valid_double(nH, 1);

    double nH0 = 0.;
    double nHp = 0.;
    double nHe0 = 0.;
    double nHep = 0.;
    double nHepp = 0.;
    double ne = 0.;

    ionization_equilibrium_calculate_densities(T, nH, X, &nH0, &nHp, &nHe0,
                                               &nHep, &nHepp, &ne);
    check_valid_double(nH0, 1);
    check_valid_double(nHp, 1);
    check_valid_double(nHe0, 1);
    check_valid_double(nHep, 1);
    check_valid_double(nHepp, 1);
    check_valid_double(ne, 1);

    /* Check for heating rates divided by number density, which is a parameter
     * that grackle requires. Check hydrogen only, and only 1 radiation group.
     * Then try with "fully neutral" and "fully ionized" gas,
     * where we assume that "fully ionized" is taken to mean that the number
     * density of neutral hydrogen is multiplied by a factor of 1e-9.
     * The cross sections assume a blackbody spectrum of T = 1e5K. */

    const double cse = 1.096971e-18;
    const double csn = 1.630511e-18;
    const double mean_photon_energy = 4.744e-11; /* erg */
    const double E_ion = 2.179e-11; /* erg, ionization energy of hydrogen */

    const double Eic = radiation_energy_density * const_speed_light_c;
    const double Nic = Eic / mean_photon_energy;

    check_valid_double(Eic, 0);
    check_valid_double(Nic, 0);

    const double number_density_units = units->density_units / const_mh;
    check_valid_double(number_density_units, 0);
    const double nHI_cgs = nH0 * number_density_units;
    check_valid_double(nHI_cgs, 1);

    const double HI_heating_rate =
        (cse * mean_photon_energy - E_ion * csn) * Nic * nHI_cgs;
    check_valid_double(HI_heating_rate, 0);
    const double HI_ionization_rate = csn * Nic;
    check_valid_double(HI_ionization_rate, 0);

    /* heating rate is expected to be in erg/s/cm^3 / nHI*/
    const double HI_heating_rate_for_grackle = HI_heating_rate / nHI_cgs;
    check_valid_double(HI_heating_rate_for_grackle, 1);
    /* ion. rate is expected to be in 1/time_units, so divide by (1/time_units)
     */
    const double HI_ionization_rate_for_grackle =
        HI_ionization_rate * units->time_units;
    check_valid_double(HI_ionization_rate_for_grackle, 1);

    if (verbose) {
      message("Heating rate: %.4e Ionization rate: %.4e", HI_heating_rate,
              HI_ionization_rate);
      message("Heating rate GR: %.4e Ionization rate GR: %.4e",
              HI_heating_rate_for_grackle, HI_ionization_rate_for_grackle);
    }
    heating_rates[x] = HI_heating_rate_for_grackle;

    /* Now check that the actual computation within grackle doesn't exceed
     * limits */
    /* using volumetric heating rate */
    double edot_volumetric = HI_heating_rate / coolunit / (dom * dom);
    check_valid_double(edot_volumetric, 1);
    double HI = nHI_cgs * const_mh / units->density_units;
    check_valid_double(HI, 1);
    double edot_RT_heating = HI_heating_rate_for_grackle / coolunit * HI / dom;
    check_valid_double(edot_RT_heating, 1);
    if (verbose)
      message("Edot volumetric: %.6e RT: %.6e", edot_volumetric,
              edot_RT_heating);
  }

  check_doubles_equal(heating_rates[0], heating_rates[1]);
  check_doubles_equal(heating_rates[1], heating_rates[2]);
}

/**
 * @brief Check whether other gas quantities derived from the
 * initial conditions are within an acceptable range
 *
 * @param radEnergy radiation energy to use
 * @param radName name of the test case. NO SPACES.
 * @param density gas density to use
 * @param name name of the test case. NO SPACES.
 * @param T temperature to deal with
 * @param params simulation parameters
 * @param units the internal units to be used in the run
 * @param verbose are we talkative?
 **/
void check_radiation_energies(float radEnergy, char *radName, float density,
                              char *name, float T, const struct simulation_params *params, const struct units* units, int verbose) {

  char fullname[80];
  const double mean_partV = params->boxsize * params->boxsize * params->boxsize / (double)params->npart;
  double volumes[3] = {mean_partV, 1.e-3 * mean_partV, 1.e3 * mean_partV};
  char *vnames[3] = {"V_av", "V_min", "V_max"};

  if (radEnergy == 0.) {
    sprintf(fullname, "E_rad:%s, %s", radName, name);
    message("radiation energy = 0 case %s; skipping.", fullname);
    return;
  }

  for (int v = 0; v < 3; v++) {
    sprintf(fullname, "E_rad:%s, %s, %s", radName, name, vnames[v]);
    const double partV = volumes[v];
    const double rad_energy_density = radEnergy / partV;
    check_valid_float(rad_energy_density, 0);

    check_grackle_internals(density, rad_energy_density, fullname, T, units, verbose);
  }
}

/**
 * @brief Check whether other gas quantities derived from the
 * initial conditions are within an acceptable range
 *
 * @param density gas density to use
 * @param name name of the test case. NO SPACES.
 * @param T temperature to deal with
 * @param units internal units to run simulation with
 * @param verbose are we talkative?
 **/
void check_luminosities(float luminosity, float density, char *name, float T,
                        const struct simulation_params *p, const struct units* units,
                        int verbose) {

  char fullname[80];
  sprintf(fullname, "luminosities: L=%.3g, %s", luminosity, name);

  if (luminosity == 0.) {
    message("luminosity = 0 case %s; skipping.", fullname);
    return;
  }

  float rad_energy_density = conversions_radiation_energy_density_from_luminosity(luminosity, p, units);
  check_grackle_internals(density, rad_energy_density, fullname, T, units, verbose);
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
  /* This needs to be the parameter filename that you plan
   * on running SWIFT with for your simulation */
  /* char *sim_run_params_filename = "swift_parameters.yml"; */
  /* char *sim_run_params_filename = "ilievTest0part2.yml"; */
  char *sim_run_params_filename = "rt_advection1D_high_redshift.yml";
  /* This is the parameter filename for the params that you
   * either set up manually or extracted from the ICs using
   * the provided script. */
  char *IC_params_filename = "advectionProblemICs.yml";

  /* Print a lot of information to the screen? */
  int verbose = 0;

  /* Initialize all sorts of variables and structs. */
  struct simulation_params simulation_params;
  simulation_params_init(&simulation_params);
  struct units units;
  units_init(&units);

  /* Read the SWIFT parameter file, and the parameters */
  /* Note that this is the original struct from swift. This way,
   * we can use the parser that came along with it. */
  struct swift_params *swift_params =
      (struct swift_params *)malloc(sizeof(struct swift_params));

  params_read_paramfile(swift_params, sim_run_params_filename);
  params_read_swift_params(swift_params, &units, &simulation_params);

  /* Get additional internal unit conversions before we continue */
  units_get_internal_units(&units);

  /* Read the simulation data parameter file, and the parameters */
  struct swift_params *sim_params =
      (struct swift_params *)malloc(sizeof(struct swift_params));
  params_read_paramfile(sim_params, IC_params_filename);
  params_read_ic_params(sim_params, &units, &simulation_params, verbose);

  /* Print out the units and parameters for a visual inspection */
  simulation_params_print(&units, &simulation_params);


  /* ------------------------*/
  /* Prepare to run examples */
  /* ------------------------*/

  /* If the min/max densities are too close to the average, re-size them */
  if (fabs(1. -simulation_params.density_min /simulation_params.density_average) < 0.05) {
    message("density_min too close to average. Resizing %.3e -> %.3e",
          simulation_params.density_min, 0.8 *simulation_params.density_average);
    simulation_params.density_min = 0.8 *simulation_params.density_average;
    check_valid_float(simulation_params.density_min, 1);
  }
  if (fabs(simulation_params.density_max / simulation_params.density_average - 1.) < 0.05) {
    message("density_max too close to average. Resizing %.3e -> %.3e",
            simulation_params.density_max, 1.2 * simulation_params.density_average);
    simulation_params.density_max = 1.2 * simulation_params.density_average;
    check_valid_float(simulation_params.density_max, 1);
  }

  /* If the min/max densities are too close to the average, re-size them */
  if (simulation_params.rad_energy_av > 0.) {
    if (simulation_params.rad_energy_min / simulation_params.rad_energy_av > 1.e-3) {
      message("photon energy_min too close to average. Resizing %.3e -> %.3e",
              simulation_params.rad_energy_min, 1.e-3 * simulation_params.rad_energy_av);
      simulation_params.rad_energy_min = 1.e-3 * simulation_params.rad_energy_av;
      check_valid_float(simulation_params.rad_energy_min, 0);
    }
    if (simulation_params.rad_energy_max / simulation_params.rad_energy_av < 1.e3) {
      message("photon energy_max too close to average. Resizing %.3e -> %.3e",
              simulation_params.rad_energy_max, 1.e3 * simulation_params.rad_energy_av);
      simulation_params.rad_energy_max = 1.e3 * simulation_params.rad_energy_av;
      check_valid_float(simulation_params.rad_energy_max, 0);
    }
  }

  /* Set up Initial conditions for radiation */
#if RT_NGROUPS == 3
  /* Fixed luminosity to heat the gas. In erg / cm^2 / s */
  double L_heating_test_cgs[3] = {1.350e+01, 2.779e+01, 6.152e+00};
#elif RT_NGROUPS == 1
  /* Fixed luminosity to heat the gas. In erg / cm^2 / s */
  double L_heating_test_cgs[1] = {4.774e+01};
#else
#error Only RT_NGROUPS = [1, 3] implemented for now
#endif
  double Erad_heating_test_cgs[RT_NGROUPS];
  for (int g = 0; g < RT_NGROUPS; g++) {
    Erad_heating_test_cgs[g] = L_heating_test_cgs[g] / const_speed_light_c;
    check_valid_double(Erad_heating_test_cgs[g], 0);
  }

  /* energy densities from stellar luminosities */
  double Erad_luminosity_test_cgs[RT_NGROUPS];
  for (int g = 0; g < RT_NGROUPS; g++) {
    Erad_luminosity_test_cgs[g] =
        conversions_radiation_energy_density_from_luminosity(simulation_params.star_emission_rates[g], &simulation_params,  &units);
    Erad_luminosity_test_cgs[g] *= units.energy_density_units;
    check_valid_double(Erad_luminosity_test_cgs[g], 0);
  }
  /* Set up arrays to loop over */
  float dens_arr[3] = {simulation_params.density_average, simulation_params.density_min, simulation_params.density_max};
  char *dens_names[3] = {"rho_av", "rho_min", "rho_max"};
  float T_test[7] = {10., 100., 1000., 1.e4, 1.e5, 1.e6, 1.e7};

  float rad_arr[3] = {simulation_params.rad_energy_av, simulation_params.rad_energy_min, simulation_params.rad_energy_max};
  char *rad_names[3] = {"Erad_av", "Erad_min", "Erad_max"};

  /* Run Gas Quantities Checks */
  /* ------------------------- */
  /* Set up some radiation energy density. 4.774e+01 erg/s/cm^2 corresponds
   * to a blackbody spectrum with T=10^5K and Ndot = 10^12 photons/s */
  const double fixed_luminosity_gas_check_cgs = 4.774e+01; /* erg/s/cm^2 */
  const double Ei = fixed_luminosity_gas_check_cgs / const_speed_light_c;
  for (int d = 0; d < 3; d++) {
    float rho = dens_arr[d];
    char *name = dens_names[d];
    for (int t = 0; t < 7; t++) {
      float T = T_test[t];
      check_gas_quantities(rho, name, T, &simulation_params, &units, verbose);
      check_grackle_internals(rho, Ei, name, T, &units, verbose);
    }
  }

  /* Run Radiation Quantities Checks */
  /* ------------------------------- */
  /* Check radiation energies only if there are some present
   * in the ICs. The radiation energies resulting from injection
   * from stars will be checked with check_luminosities() */
  for (int d = 0; d < 3; d++) {
    float rho = dens_arr[d];
    char *name = dens_names[d];

    for (int t = 0; t < 7; t++) {
      float T = T_test[t];

      for (int e = 0; e < 3; e++) {
        float radEnergy = rad_arr[e];
        char *radName = rad_names[e];
        if (simulation_params.rad_energy_av > 0.)
          check_radiation_energies(radEnergy, radName, rho, name, T, &simulation_params, &units, verbose);
      }

      for (int g = 0; g < RT_NGROUPS; g++) {
        float l = simulation_params.star_emission_rates[g];
        check_luminosities(l, rho, name, T, &simulation_params, &units, verbose);
      }
    }
  }

  /* Run Grackle cooling test */
  /* ------------------------ */
  for (int d = 0; d < 3; d++) {
    float rho = dens_arr[d];
    char *name = dens_names[d];
    run_grackle_cooling_test(rho, name, units.mass_units, units.length_units, units.time_units, units.density_units, units.velocity_units, units.internal_energy_units, verbose);
  }

  /* Run Grackle heating test */
  /* ------------------------ */
  for (int d = 0; d < 3; d++) {
    float rho = dens_arr[d];
    char *name = dens_names[d];
    run_grackle_heating_test(rho, Erad_heating_test_cgs, name, units.mass_units,
                             units.length_units, units.time_units, units.density_units,
                             units.velocity_units, units.internal_energy_units,
                             /*dump_results=*/1, verbose);
  }

  /* Run Grackle heating test with provided luminosities */
  /* --------------------------------------------------- */

  if (simulation_params.use_const_emission_rates) {
    for (int d = 0; d < 3; d++) {
      float rho = dens_arr[d];
      char fullname[80];
      sprintf(fullname, "%s-NO_OUTPUT-LUMINOSITY_TEST", dens_names[d]);
      run_grackle_heating_test(
          rho, Erad_luminosity_test_cgs, fullname, units.mass_units, units.length_units,
          units.time_units, units.density_units, units.velocity_units, units.internal_energy_units,
          /*dump_results=*/0, verbose);
    }
  } else {
    error("This test is not set up without constant emission rates");
  }

  /* Clean up after yourself */
  free(swift_params);
  free(sim_params);

  message("Check completed with %d warning(s).", warnings);
  message("Bye, and good luck.");

  return 0;
}
