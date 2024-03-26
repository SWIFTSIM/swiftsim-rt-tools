/* ---------------------------------------------
 * In this example, we start with high internal
 * energies and a fully ionized gas, and just
 * let it cool without any RT.
 * --------------------------------------------- */

/* define this before including my_grackle_utils.h */
#define FIELD_SIZE 1
/* Skip grackle warnings of deprecated local functions. */
#define OMIT_LEGACY_INTERNAL_GRACKLE_FUNC

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <grackle.h>

#include "constants.h"
#include "cosmology.h"
#include "ionization_equilibrium.h"
#include "mean_molecular_weight.h"
#include "my_grackle_utils.h"

int main() {

  /******************************************************************
   * Set up initial conditions and runtime parameters.
   *****************************************************************/

  /* Print some extra data to screen? */
  int verbose = 1;
  /* output frequency in number of steps */
  const int output_frequency = 1;
  /* Integrate in intervals of dlog a ? */
  const int log_integration = 0;
  /* How many steps to run */
  const int nsteps = 20000;
  /* Option to turn off cosmological integration. Intended to comare results
   * of cosmo and non-cosmo outputs with identical ICs run over the identical
   * time frame, specified by a_begin and a_end.
   * Ideally, we want to find an example where the results are visibly
   * different. */
  const int with_cosmo = 1;
  /* output file */
  FILE *fd;
  if (with_cosmo) {
    fd = fopen("out.dat", "w");
  } else {
    fd = fopen("outNoCosmo.dat", "w");
  }

  /* Define units : use the same as internal units for swift */
  /* ------------------------------------------------------- */
  double mass_units = 1.99848e43;
  double length_units = 3.08567758e21;
  double velocity_units = 1e5;

  double density_units =
      mass_units / length_units / length_units / length_units;
  double time_units = length_units / velocity_units;

  /* Cosmology                  */
  /* -------------------------- */
  /* double a_begin = 0.0099;  [> z~100 <] */
  /* double a_end = 0.014081;  [> z~70 <] */
  /* double a_begin = 0.0625;  [> z~15 <] */
  /* double a_end = 0.09091;  [> z~10 <] */
  double a_begin = 0.0476; /* z~20 */
  /* double a_begin = 0.09091;  [> z~10 <] */
  /* double a_end = 0.166667; [> z~5 <] */
  double a_end = 0.2; /* z~4 */
  /* double a_end = 1.0; [> z~5 <] */

  const double log_a_begin = log(a_begin);
  const double log_a_end = log(a_end);
  /* Only one of these will be used, depending on whether you set
   * int log_integration = 1 */
  const double dlog_a = (log_a_end - log_a_begin) / nsteps;
  const double da = (a_end - a_begin) / nsteps;

  /* Use this a for conversions to comoving frame. */
  double a_convert_comoving = 1.;
  if (with_cosmo)
    a_convert_comoving = a_begin;

  struct cosmology cosmology;

  /* Planck13 (EAGLE flavour) */
  cosmology.Omega_cdm = 0.2587; /* Dark matter density parameter*/
  cosmology.Omega_b = 0.04825;  /* baryon density parameter*/
  cosmology.Omega_l = 0.693;    /* Dark Energy density parameter */
  cosmology.Omega_k = 0.;       /* Radiation density parameter */
  cosmology.Omega_r = 0.;       /* Radiation density parameter */
  cosmology.Omega_nu = 0.;      /* Neutrino density parameter */
  cosmology.w_0 = -1.0; /* Dark-energy equation-of-state parameter at z=0. */
  cosmology.w_a =
      0.; /* Dark-energy equation-of-state time evolution parameter. */
  cosmology.H_0 = 67.77; /* Hubble constant at z=0 in km/s/Mpc */

  cosmo_convert_H0_to_internal_units(&cosmology, time_units);

  /* Compute a(t) and t(a) tables for interpolation */
  /* ---------------------------------------------- */

  double a_table[COSMO_TABLE_ELEMENTS];
  double t_table[COSMO_TABLE_ELEMENTS];

  cosmo_get_tables(a_table, t_table, &cosmology, a_begin, a_end);

  /* Set up initial conditions for gas cells */
  /* --------------------------------------- */
  double hydrogen_fraction_by_mass = 0.76;
  /* Use solution of swift's output. This is in internal units already. */
  /* However, these are in physical units. We'll turn them into comoving units
   * later. */
  double gas_density_phys = 0.00024633363 / 20.;
  double internal_energy_phys = 21201.879;

  /* Derived quantities from ICs */
  /* --------------------------- */

  /* Assuming fully ionized gas */
  double mu_init = mean_molecular_weight_from_mass_fractions(
      0., hydrogen_fraction_by_mass, 0., 0., (1. - hydrogen_fraction_by_mass));
  double internal_energy_phys_cgs = internal_energy_phys * length_units *
                                    length_units / (time_units * time_units);

  double T_phys = internal_energy_phys_cgs * (const_adiabatic_index - 1) *
                  mu_init * const_mh / const_kboltz;
  if (verbose)
    printf("Initial setup: u_cgs %g T_cgs %g [physical units]\n",
           internal_energy_phys_cgs, T_phys);

  /* define the hydrogen number density */
  /* use `gr_float` to use the same type of floats that grackle
   * is compiled in. Compile grackle with precision-32 if you want
   * floats, or precision-64 otherwise. */
  gr_float nH_phys =
      hydrogen_fraction_by_mass * gas_density_phys / (const_mh / mass_units);
  gr_float nHI;
  gr_float nHII;
  gr_float nHeI;
  gr_float nHeII;
  gr_float nHeIII;
  gr_float ne;

  /* get densities of primordial species assuming ionization equilibrium */
  ionization_equilibrium_calculate_densities(
      T_phys, nH_phys, hydrogen_fraction_by_mass, &nHI, &nHII, &nHeI, &nHeII,
      &nHeIII, &ne);

  gr_float HI_density = nHI * (const_mh / mass_units);
  gr_float HII_density = nHII * (const_mh / mass_units);
  gr_float HeI_density = nHeI * (4 * const_mh / mass_units);
  gr_float HeII_density = nHeII * (4 * const_mh / mass_units);
  gr_float HeIII_density = nHeIII * (4 * const_mh / mass_units);
  /* !! this is the convention adopted by Grackle
   * !! e_density is the electron number density multiplied by proton mass,
   * !! or electron mass density * nH / ne */
  gr_float e_density = ne * (const_mh / mass_units);

  /* Convert them to co-moving frame*/
  HI_density = cosmo_get_comoving_density(HI_density, a_convert_comoving);
  HII_density = cosmo_get_comoving_density(HII_density, a_convert_comoving);
  HeI_density = cosmo_get_comoving_density(HeI_density, a_convert_comoving);
  HeII_density = cosmo_get_comoving_density(HeII_density, a_convert_comoving);
  HeIII_density = cosmo_get_comoving_density(HeIII_density, a_convert_comoving);
  e_density = cosmo_get_comoving_density(e_density, a_convert_comoving);

  /* Store them all in a single array for simplicity. */
  gr_float species_densities[12] = {
      HI_density, HII_density, HeI_density, HeII_density, HeIII_density,
      e_density,  0.,          0.,          0.,           0.,
      0.,         0.};

  /* Grackle behaviour setup */
  /* ----------------------- */

  int UVbackground = 0;         /* toogle on/off the UV background */
  int primordial_chemistry = 1; /* choose the chemical network */
  int use_radiative_cooling = 1;
  int use_radiative_transfer = 0;
  char *grackle_data_file = "";

  /* These are in cgs anyway, so no transform to co-moving coordinates
   * necessary. */
  gr_float RT_HI_ionization_rate = 0.;
  gr_float RT_HeI_ionization_rate = 0.;
  gr_float RT_HeII_ionization_rate = 0.;
  gr_float RT_H2_dissociation_rate = 0.;
  gr_float RT_heating_rate = 0.;

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
  setup_grackle_units_cosmo(&grackle_units_data, density_units, length_units,
                            time_units, a_begin, with_cosmo);

  /* Chemistry Parameters */
  /* -------------------- */

  /* Second, create a chemistry object for parameters.  This needs to be a
   * pointer. */
  chemistry_data grackle_chemistry_data;
  if (local_initialize_chemistry_parameters(&grackle_chemistry_data) == 0) {
    fprintf(stderr, "Error in local_initialize_chemistry_parameters.\n");
    return EXIT_FAILURE;
  }

  /* Set parameter values for chemistry. */
  setup_grackle_chemistry(&grackle_chemistry_data, primordial_chemistry,
                          UVbackground, grackle_data_file,
                          use_radiative_cooling, use_radiative_transfer,
                          hydrogen_fraction_by_mass);

  /* Initialize the chemistry_data_storage object to be able to use local
   * functions */
  chemistry_data_storage grackle_chemistry_rates;
  if (local_initialize_chemistry_data(&grackle_chemistry_data,
                                      &grackle_chemistry_rates,
                                      &grackle_units_data) == 0) {
    fprintf(stderr, "Error in local_initialize_chemistry_data.\n");
    return EXIT_FAILURE;
  }

  /* Gas Data */
  /* -------- */

  /* Create struct for storing grackle field data */
  double gas_density_comoving =
      cosmo_get_comoving_density(gas_density_phys, a_convert_comoving);
  double internal_energy_comoving = cosmo_get_comoving_internal_energy(
      internal_energy_phys, a_convert_comoving);

  grackle_field_data grackle_fields;
  setup_grackle_fields(&grackle_fields, species_densities, interaction_rates,
                       gas_density_comoving, internal_energy_comoving);

  /* Write headers */
  /* ------------- */

  /* First to stdout */

  if (verbose) {
    printf("%15s%15s%15s%15s%15s%15s%15s%15s\n",
           "Initial setup: ", "Temperature", "nHI", "nHII", "nHeI", "nHeII",
           "nHeIII", "ne");
    printf("%15s%15g%15g%15g%15g%15g%15g%15g\n\n", "Initial setup: ", T_phys,
           nHI, nHII, nHeI, nHeII, nHeIII, ne);
  }

  write_cosmo_header(stdout);
  write_cosmo_timestep(stdout, &grackle_fields, &grackle_units_data,
                       &grackle_chemistry_data, &grackle_chemistry_rates,
                       /*field_index=*/0, /*t=*/0., /*dt=*/0., a_begin,
                       time_units, /*step=*/0, with_cosmo);

  /* And now to the output file. */
  write_my_cosmo_setup(fd, grackle_fields, &grackle_chemistry_data, mass_units,
                       length_units, velocity_units, a_begin, a_end, &cosmology,
                       hydrogen_fraction_by_mass, gas_density_phys,
                       internal_energy_phys, with_cosmo);
  write_cosmo_header(fd);
  write_cosmo_timestep(fd, &grackle_fields, &grackle_units_data,
                       &grackle_chemistry_data, &grackle_chemistry_rates,
                       /*field_index=*/0, /*t=*/0., /*dt=*/0., a_begin,
                       time_units, /*step=*/0, with_cosmo);

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/

  int step = 0;
  /* Current expansion scale factor */
  double a = a_begin;
  /* log of current expansion scale factor */
  double log_a = log_a_begin;
  /* Current time */
  double t = 0.;
  /* Current time step size */
  double dt = 0.;

  /* Value of expansion factor at the end of the step */
  double a_next = a_begin;

  while (a < a_end) {

    a = a_next;
    t += dt;
    step += 1;

    /* Compute next time step size. */
    if (log_integration) {
      /* Marching in steps of equal dlog(a) */
      log_a += dlog_a;
      a_next = exp(log_a);
    } else {
      /* Marching in steps of equal a */
      a_next += da;
    }

    /* Apply change in a to grackle. */
    update_grackle_units_cosmo(&grackle_units_data, density_units, length_units,
                               a, with_cosmo);

    /* Get time step size. */
    dt = cosmo_get_dt(a, a_next, a_begin, a_end, t_table);

    /* Solve chemistry. */
    if (local_solve_chemistry(&grackle_chemistry_data, &grackle_chemistry_rates,
                              &grackle_units_data, &grackle_fields, dt) == 0) {
      fprintf(stderr, "Error in solve_chemistry.\n");
      return EXIT_FAILURE;
    }

    write_cosmo_timestep(stdout, &grackle_fields, &grackle_units_data,
                         &grackle_chemistry_data, &grackle_chemistry_rates,
                         /*field_index=*/0, t, dt, a, time_units, step,
                         with_cosmo);

    if (step % output_frequency == 0)
      write_cosmo_timestep(fd, &grackle_fields, &grackle_units_data,
                           &grackle_chemistry_data, &grackle_chemistry_rates,
                           /*field_index=*/0, t, dt, a, time_units, step,
                           with_cosmo);
  }

  /* Clean up after yourself */

  fflush(stdout);
  fclose(fd);
  clean_up_fields(&grackle_fields);
  local_free_chemistry_data(&grackle_chemistry_data, &grackle_chemistry_rates);

  return EXIT_SUCCESS;
}
