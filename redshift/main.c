/* ---------------------------------------------
 * In this example, we start with high internal
 * energies and a fully ionized gas, and just
 * let it cool without any RT.
 * --------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "cosmology.h"
#include "redshift.h"

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
  const int nsteps = 1000;
  /* Option to turn off cosmological integration. Intended to comare results
   * of cosmo and non-cosmo outputs with identical ICs run over the identical
   * time frame, specified by a_begin and a_end.
   * Ideally, we want to find an example where the results are visibly
   * different. */
  const int with_cosmo = 1;
  /* Option to turn off shifting of the BB spectrum using T_bb ~ 1/a */
  const int with_shifting_BB = 1;
  const int flat_spectrum = 0;
  /* output file */
  FILE *fd;
  if (with_shifting_BB) {
    fd = fopen("out.dat", "w");
  } else {
    fd = fopen("out.dat", "w");
  }

  /* Define units : use the same as internal units for swift */
  /* ------------------------------------------------------- */
  double mass_units = 1.99848e43;
  double length_units = 3.08567758e21;
  double velocity_units = 1e5;

  double density_units =
      mass_units / length_units / length_units / length_units;
  double time_units = length_units / velocity_units;
  double energy_units = mass_units * velocity_units * velocity_units;
  (void)energy_units;
  (void)density_units;

  /* Cosmology                  */
  /* -------------------------- */
  /* double a_begin = 0.0099;  [> z~100 <] */
  /* double a_end = 0.014081;  [> z~70 <] */
  /* double a_begin = 0.0625;  [> z~15 <] */
  /* double a_end = 0.09091;  [> z~10 <] */
  /* double a_begin = 0.0476; [> z~20 <] */
  /* double a_end = 0.166667; [> z~5 <] */
  double a_begin = 0.0909;
  double a_end = 0.1428571;

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
  cosmology.Omega_cdm = 0.2587481; /* Dark matter density parameter*/
  cosmology.Omega_b = 0.0482519;  /* baryon density parameter*/
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
  /* Use solution of swift's output. This is in internal units already. */
  /* However, these are in physical units. We'll turn them into comoving units
   * later. */
  double gas_density_phys = 2.463186927340698e-04;
  double internal_energy_phys = 2.120106942539497e+04;
  double gas_volume_phys = 1.e10; /* Internal units */

  /* Derived quantities from ICs */
  /* --------------------------- */

  /* Gas Data */
  /* -------- */

  /* Read in photon energy from file */
  FILE *energy_file = fopen("photon_energies", "rb");

  double photon_energy[RT_NGROUPS];
  size_t temp = fread(photon_energy, sizeof(double), LEN(photon_energy), energy_file);
  (void)temp;
  fclose(energy_file); 

  for (int g = 0; g < RT_NGROUPS; g++) {
    if (flat_spectrum) {
      photon_energy[g] = 1;
    } else {
      photon_energy[g] *= 1e-19;
    }
    //printf("%.10e\n", photon_energy[g]);
  }
  /* Read photon frequency bins from file */
  FILE *bins_file = fopen("frequency_bins", "rb");

  double frequency_bins[RT_NGROUPS];
  temp = fread(frequency_bins, sizeof(double), LEN(frequency_bins), bins_file);
  (void)temp;
  fclose(bins_file);

  double BB_temperature_phys = 1e6; /* Kelvin */
  // (void)BB_temperature;
  // (void)frequency_bins;


  /* Create struct for storing grackle field data */
  double gas_density_comoving =
      cosmo_get_comoving_density(gas_density_phys, a_convert_comoving);
  double internal_energy_comoving = cosmo_get_comoving_internal_energy(
      internal_energy_phys, a_convert_comoving);
  double gas_volume_comoving = cosmo_get_comoving_volume(gas_volume_phys, a_convert_comoving);
  double BB_temperature_comoving = cosmo_get_comoving_radiative_temperature(BB_temperature_phys, a_convert_comoving);
  (void)gas_density_comoving;
  (void)internal_energy_comoving;
  
  /* Write headers */
  /* ------------- */

  /* First to stdout */

  double total_photon_energy = 0.;
  for (int i = 0; i < RT_NGROUPS; i++) {
    total_photon_energy += photon_energy[i];
  }
  
  if (verbose) {
    printf("%15s%15s\n",
           "Initial setup: ", "Photon energy");
    printf("%15s%15f\n\n", "Initial setup: ", total_photon_energy);
  }

  write_cosmo_header(stdout);
  write_cosmo_timestep(stdout, 0, a_begin, total_photon_energy, gas_volume_phys, BB_temperature_phys, with_shifting_BB);

  write_cosmo_header(fd);
  write_cosmo_timestep(fd, 0, a_begin, total_photon_energy, gas_volume_phys, BB_temperature_phys, with_shifting_BB);

  /*********************************************************************
  / Calling the redshift solver
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


  double redshift_factor[RT_NGROUPS];

  while (a < a_end) {    
    
    /* Reset total photon energy */
    total_photon_energy = 0.;

    /* Recalculate BB temperature */
    double new_BB_temperature;
    if (with_shifting_BB) {
      new_BB_temperature = cosmo_get_physical_radiative_temperature(BB_temperature_comoving, a);
    } else {
      new_BB_temperature = BB_temperature_phys;
    }
    calculate_redshift_factor(redshift_factor, frequency_bins, new_BB_temperature);
    
    
    double H = 1 / (time_integrand(a, &cosmology) * a);

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
    
    dt = cosmo_get_dt(a, a_next, a_begin, a_end, t_table);

    /* Update photon energy */
    for (int i = 0; i < RT_NGROUPS; i++) {
      printf("[group %d] - Redshift factor : %.10e\n", i+1, redshift_factor[i]);
      //printf("[group %d] - Photon energy   : %.10e\n", i+1, photon_energy[i]);
      /*if (redshift_factor[i] > 0) {
	if (i < RT_NGROUPS - 1) {
          photon_energy[i] -= H * dt * (redshift_factor[i] * photon_energy[i+1] + 0*photon_energy[i]);
	}
      } else {
	photon_energy[i] += H * dt * (redshift_factor[i] * photon_energy[i] + 0*photon_energy[i]);
      }
      */ 
      
      photon_energy[i] -= H * dt * (1 * photon_energy[i] + 0*photon_energy[i]);
      /*if (i+1 < RT_NGROUPS) {
        photon_energy[i] += H * dt * (redshift_factor[i] * photon_energy[i+1] + 0*photon_energy[i+1]);
      } */
    }
    
    for (int i = 0; i < RT_NGROUPS; i++) {
      total_photon_energy += photon_energy[i];
    }

    /* Update volume */
    gas_volume_phys = cosmo_get_physical_volume(gas_volume_comoving, a);
    
    write_cosmo_timestep(stdout, step, a, total_photon_energy, gas_volume_phys, new_BB_temperature, with_shifting_BB);

    if (step % output_frequency == 0)
      write_cosmo_timestep(fd, step, a, total_photon_energy, gas_volume_phys, new_BB_temperature, with_shifting_BB);
  }

  /* Clean up after yourself */

  fflush(stdout);
  fclose(fd);

  return EXIT_SUCCESS;
}
