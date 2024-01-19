#include "conversions.h"
#include "constants.h"
#include "error.h"
#include "validity_check_macros.h"

extern int warnings;

/**
 * @brief estimate a time step size.
 * @return time step in internal units.
 **/
float conversions_estimate_dt(const struct simulation_params *p) {
  if (p->smoothing_length == 0.)
    error("sml=0?");
  if (p->c_reduced == 0.)
    error("c_red=0?");
  return p->smoothing_length / p->c_reduced;
}

/**
 * @brief Get an estimate for the injected radiation energy density
 * from a given luminosity.
 *
 * @param luminosity luminosity in units of solar luminosities
 * @return the radiation energy density in internal units.
 *
 **/
float conversions_radiation_energy_density_from_luminosity(
    double luminosity, const struct simulation_params *p,
    const struct units *units) {

  /* To obtain some energy density, set E = L * dt / V
   * If you want to vary particle volume, you should also vary dt
   * by the same factor according to the CFL logic. This means
   * however that the factor to modify dt and V cancels out since
   * we need dt/V, so no reason to check several scenarios there. */

  const double partV = p->boxsize * p->boxsize * p->boxsize / (double)p->npart;
  check_valid_float(partV, 0);
  const double dt = conversions_estimate_dt(p);
  check_valid_float(dt, 0);
  float rad_energy_density =
      luminosity * const_L_Sun / units->power_units * dt / partV;
  check_valid_float(rad_energy_density, 0);

  return rad_energy_density;
}
