#ifndef GRACKLE_CHECKS_H
#define GRACKLE_CHECKS_H

/* A variety of checks for grackle's consistency. */

#include "constants.h"
#include "error.h"
#include <grackle.h>

/**
 * @brief check that the sum of the components adds up
 * to the total density.
 * This assumes that the gas density remains constant each step.
 *
 * @param density the initial density used in test
 * @param grackle_fields the struct containing the gas fields
 */
void grackle_checks_density_sum(float density,
                                grackle_field_data *grackle_fields) {

  double dens_sum =
      grackle_fields->HI_density[0] + grackle_fields->HII_density[0] +
      grackle_fields->HeI_density[0] + grackle_fields->HeII_density[0] +
      grackle_fields->HeIII_density[0] +
      /* grackle_fields->e_density[0] + */
      grackle_fields->HM_density[0] + grackle_fields->H2I_density[0] +
      grackle_fields->H2II_density[0] + grackle_fields->DI_density[0] +
      grackle_fields->DII_density[0] + grackle_fields->HDI_density[0];

  if (fabs(grackle_fields->density[0] / dens_sum - 1.) > 0.0001)
    error("Error: Density sum does not correspond to total density; "
          "total=%.6e sum=%.6e ratio=%.3e",
          grackle_fields->density[0], dens_sum,
          grackle_fields->density[0] / dens_sum);

  if (fabsf((float)grackle_fields->density[0] / density - 1.f) > 0.0001)
    error("Error: initial and current density not equal; "
          "initial=%.6e current=%.6e ratio=%.3e",
          density, grackle_fields->density[0],
          grackle_fields->density[0] / density);
}

/**
 * @brief check that the sum of the ion numbers adds up
 * to the electron number.
 * This assumes that the gas density remains constant each step.
 *
 * @param grackle_fields the struct containing the gas fields
 */
void grackle_checks_ion_sum(grackle_field_data *grackle_fields,
                            double mass_units) {

  double mh = const_mh / mass_units;

  double ion_number_sum = grackle_fields->HII_density[0] / mh +
                          grackle_fields->HeII_density[0] / (4. * mh) +
                          2. * grackle_fields->HeIII_density[0] / (4. * mh);
  /* grackle_fields->H2II_density[0] / (2. * mh)  */
  /* grackle_fields->DII_density[0] / (2. * mh); */

  /* e_density is the electron number density multiplied by proton mass */
  double ne = grackle_fields->e_density[0] / mh;

  if (ne > TINY_NUMBER && fabs(ion_number_sum / ne - 1.) > 0.001)
    error("Error: Ion number sum does not correspond to total density; "
          "sum=%.6e ne=%.6e ratio=%.3e",
          ion_number_sum, ne, ion_number_sum / ne);
}
#endif
