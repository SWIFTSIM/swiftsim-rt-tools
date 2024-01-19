#include "units.h"
#include "validity_check_macros.h"

extern int warnings;

/**
 * Initialize the units struct.
 */
void units_init(struct units *units) {

  units->mass_units = 0.;
  units->time_units = 0.;
  units->length_units = 0.;
  units->density_units = 0.;
  units->velocity_units = 0.;
  units->temperature_units = 0.;
  units->internal_energy_units = 0.;
  units->energy_units = 0.;
  units->energy_density_units = 0.;
  units->power_units = 0.;
}

/**
 * @brief get additional internal unit conversions. While some units
 * are being read in from parameter files, I add other derived units
 * here for convenience.
 **/
void units_get_internal_units(struct units *u) {

  /* Convert units */
  const double volume_units =
      (u->length_units * u->length_units * u->length_units);
  check_valid_double(volume_units, 0);

  u->time_units = u->length_units / u->velocity_units;
  check_valid_float(u->time_units, 0);

  u->density_units = u->mass_units / volume_units;
  check_valid_float(u->density_units, 0);

  u->internal_energy_units = u->velocity_units * u->velocity_units;
  check_valid_float(u->internal_energy_units, 0);

  u->energy_units = u->mass_units * u->internal_energy_units;
  /* It turns out that the energy units are commonly above float exponent
   * limits, but that's not an issue because the energy itself can be very
   * high too. So in the end, it works out. */
  /* check_valid_float(energy_units, 0); */
  /* Nonetheless, print a warning for now. */
  if (u->energy_units != 0. &&
      (fabs(u->energy_units) > 1e30 || fabs(u->energy_units) < 1e-30)) {
    fprintf(stdout,
            "WARNING: %s:%s:%d: energy units has large exponent: %.6e\n",
            __FILE__, __FUNCTION__, __LINE__, u->energy_units);
    warnings++;
  }

  u->energy_density_units = u->energy_units / volume_units;
  check_valid_float(u->energy_density_units, 0);

  u->power_units = u->energy_units / u->time_units;
  check_valid_float(u->power_units, 0);
}
