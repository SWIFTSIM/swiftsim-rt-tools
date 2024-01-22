/* -------------------------------------------------
 * A suite of tests both with and without grackle
 * to check whether your choice of units may produce
 * bad results.
 * ------------------------------------------------- */

/* FPE stuff. Link with -lm */
/* #define _GNU_SOURCE */
/* #include <fenv.h> */

/* define RT_NGROUPS before including local headers like my_grackle_utils.h */
#include "rt_ngroups.h"

/* Grackle related macros */
#define FIELD_SIZE 1
/* don't raise warnings about deprecated grackle functions */
#define OMIT_LEGACY_INTERNAL_GRACKLE_FUNC
#include "my_grackle_utils.h"

#include "checks.h"
#include "conversions.h"
#include "cmdlineargs.h"
#include "grackle_cooling_test.h"
#include "grackle_heating_test.h"
#include "read_params.h"
#include "units.h"
#include "validity_check_macros.h"

int warnings = 0;


int main(int argc, char* argv[]) {

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
  char sim_run_params_filename[MAX_FILENAME_SIZE];
  /* This is the parameter filename for the params that you
   * either set up manually or extracted from the ICs using
   * the provided script. */
  char IC_params_filename[MAX_FILENAME_SIZE];
  read_cmdlineargs(sim_run_params_filename, IC_params_filename, argc, argv);

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
  if (fabs(1. - simulation_params.density_min /
                    simulation_params.density_average) < 0.05) {
    message("density_min too close to average. Resizing %.3e -> %.3e",
            simulation_params.density_min,
            0.8 * simulation_params.density_average);
    simulation_params.density_min = 0.8 * simulation_params.density_average;
    check_valid_float(simulation_params.density_min, 1);
  }
  if (fabs(simulation_params.density_max / simulation_params.density_average -
           1.) < 0.05) {
    message("density_max too close to average. Resizing %.3e -> %.3e",
            simulation_params.density_max,
            1.2 * simulation_params.density_average);
    simulation_params.density_max = 1.2 * simulation_params.density_average;
    check_valid_float(simulation_params.density_max, 1);
  }

  /* If the min/max densities are too close to the average, re-size them */
  if (simulation_params.rad_energy_av > 0.) {
    if (simulation_params.rad_energy_min / simulation_params.rad_energy_av >
        1.e-3) {
      message("photon energy_min too close to average. Resizing %.3e -> %.3e",
              simulation_params.rad_energy_min,
              1.e-3 * simulation_params.rad_energy_av);
      simulation_params.rad_energy_min =
          1.e-3 * simulation_params.rad_energy_av;
      check_valid_float(simulation_params.rad_energy_min, 0);
    }
    if (simulation_params.rad_energy_max / simulation_params.rad_energy_av <
        1.e3) {
      message("photon energy_max too close to average. Resizing %.3e -> %.3e",
              simulation_params.rad_energy_max,
              1.e3 * simulation_params.rad_energy_av);
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
        conversions_radiation_energy_density_from_luminosity(
            simulation_params.star_emission_rates[g], &simulation_params,
            &units);
    Erad_luminosity_test_cgs[g] *= units.energy_density_units;
    check_valid_double(Erad_luminosity_test_cgs[g], 0);
  }
  /* Set up arrays to loop over */
  float dens_arr[3] = {simulation_params.density_average,
                       simulation_params.density_min,
                       simulation_params.density_max};
  char *dens_names[3] = {"rho_av", "rho_min", "rho_max"};
  float T_test[7] = {10., 100., 1000., 1.e4, 1.e5, 1.e6, 1.e7};

  float rad_arr[3] = {simulation_params.rad_energy_av,
                      simulation_params.rad_energy_min,
                      simulation_params.rad_energy_max};
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
          check_radiation_energies(radEnergy, radName, rho, name, T,
                                   &simulation_params, &units, verbose);
      }

      for (int g = 0; g < RT_NGROUPS; g++) {
        float l = simulation_params.star_emission_rates[g];
        check_luminosities(l, rho, name, T, &simulation_params, &units,
                           verbose);
      }
    }
  }

  /* Run Grackle cooling test */
  /* ------------------------ */
  for (int d = 0; d < 3; d++) {
    float rho = dens_arr[d];
    char *name = dens_names[d];
    run_grackle_cooling_test(rho, name, units.mass_units, units.length_units,
                             units.time_units, units.density_units,
                             units.velocity_units, units.internal_energy_units,
                             verbose);
  }

  /* Run Grackle heating test */
  /* ------------------------ */
  for (int d = 0; d < 3; d++) {
    float rho = dens_arr[d];
    char *name = dens_names[d];
    run_grackle_heating_test(rho, Erad_heating_test_cgs, name, units.mass_units,
                             units.length_units, units.time_units,
                             units.density_units, units.velocity_units,
                             units.internal_energy_units,
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
          rho, Erad_luminosity_test_cgs, fullname, units.mass_units,
          units.length_units, units.time_units, units.density_units,
          units.velocity_units, units.internal_energy_units,
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
