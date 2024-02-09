#include "checks.h"
#include "conversions.h"
#include "constants.h"
#include "validity_check_macros.h"
#include "error.h"

#define OMIT_LEGACY_INTERNAL_GRACKLE_FUNC
#include "ionization_equilibrium.h"

extern int warnings;


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
void check_gas_quantities(float density, char *name, float T,
                          const struct simulation_params *params,
                          const struct units *units, int verbose) {

  /* assume mean molecular weight of 1 for this test. While that isn't correct,
   * it should do the trick for the purpose of this test. */
  message("Checking T=%.1e case=%s", T, name);

  const float gamma = const_adiabatic_index;
  const float gamma_minus_one = gamma - 1.f;

  /* Get and check internal energy */
  const float mu = 1.;
  const float internal_energy_cgs =
      const_kboltz * T / (gamma_minus_one * mu * const_mh);
  const float internal_energy =
      internal_energy_cgs / units->internal_energy_units;
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
 * @param params simulation parameters.
 * @param units the internal units to be used in the simulation
 * @param verbose are we talkative?
 **/
void check_grackle_internals(float density, float radiation_energy_density,
                             char *name, float T,
                             const struct simulation_params *params,
                             const struct units *units,
                             int verbose) {

  message("checking %s, T=%.1e", name, T);

  /* Check that the total number density is within the valid limits for
   * grackle to handle. */
  /* We need to check the physical values here, so convert cosmo quantities
   * to physical ones. */

  double scales_arr[2] = {params->a_begin, params->a_end};

  for (int s = 0; s < 2; s++){
    const double a = scales_arr[s];
    const double scale = 1. / (a * a * a);
    const double n_cgs = density * units->density_units * scale / const_mh;
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
  }

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
                              char *name, float T,
                              const struct simulation_params *params,
                              const struct units *units, int verbose) {

  char fullname[80];
  const double mean_partV = params->boxsize * params->boxsize *
                            params->boxsize / (double)params->npart;
  double volumes[3] = {mean_partV, 1.e-3 * mean_partV, 1.e3 * mean_partV};
  char *vnames[3] = {"V_av", "V_min", "V_max"};

  if (radEnergy == 0.) {
    sprintf(fullname, "E_rad:%s, %s", radName, name);
    message("radiation energy = 0 case %s; skipping.", fullname);
    return;
  }

  double flux_factor = params->fc * const_speed_light_c / units->velocity_units;

  for (int v = 0; v < 3; v++) {
    sprintf(fullname, "E_rad:%s, %s, %s", radName, name, vnames[v]);
    const double partV = volumes[v];
    const double rad_energy_density = radEnergy / partV;
    check_valid_float(rad_energy_density, 0);

    const double flux = flux_factor * rad_energy_density;
    check_valid_float(flux, 0);
    // We do intermediate comutations containing F^2, so that needs to fit as well
    const double f2 = flux * flux;
    check_valid_float(f2, 0);

    check_grackle_internals(density, rad_energy_density, fullname, T, params, units, verbose);
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
                        const struct simulation_params *params,
                        const struct units *units, int verbose) {

  char fullname[80];
  sprintf(fullname, "luminosities: L=%.3g, %s", luminosity, name);

  if (luminosity == 0.) {
    message("luminosity = 0 case %s; skipping.", fullname);
    return;
  }

  float rad_energy_density = conversions_radiation_energy_density_from_luminosity(luminosity, params, units);
  check_grackle_internals(density, rad_energy_density, fullname, T, params, units, verbose);
  check_radiation_energies(rad_energy_density, "luminosityTest", density, name, T, params, units, verbose);

}


