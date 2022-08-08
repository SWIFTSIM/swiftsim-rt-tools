#ifndef MEAN_MOLECULAR_WEIGHT_H
#define MEAN_MOLECULAR_WEIGHT_H

#include "constants.h"
#include <grackle.h>
#include <stdlib.h>

/*********************************************************************
/ Simple function to estimate the mean molecular weight
*********************************************************************/

/* Get mean molecular weight for given species densities */
gr_float mean_weight_from_densities(gr_float density, gr_float nHI,
                                    gr_float nHII, gr_float nHeI,
                                    gr_float nHeII, gr_float nHeIII,
                                    gr_float ne, gr_float mass_units) {

  gr_float n;
  gr_float mu;

  n = nHI + nHII + nHeI + nHeII + nHeIII + ne;
  mu = density / n / (const_mh / mass_units);

  return mu;
}

/* Get mean molecular weight for given species mass fractions */
gr_float mean_molecular_weight_from_mass_fractions(gr_float XHI, gr_float XHII,
                                                   gr_float XHeI,
                                                   gr_float XHeII,
                                                   gr_float XHeIII) {

  /* 1/mu = sum_j X_j / A_j * (1 + E_j)
   * A_H    = 1, E_H    = 0
   * A_Hp   = 1, E_Hp   = 1
   * A_He   = 4, E_He   = 0
   * A_Hep  = 4, E_Hep  = 1
   * A_Hepp = 4, E_Hepp = 2             */
  const gr_float one_over_mu =
      XHI + 2.0 * XHII + 0.25 * XHeI + 0.5 * XHeII + 0.75 * XHeIII;

  return (1.0 / one_over_mu);
}

/*********************************************************************
/ Compute the mean molecular weight using the same rules as Grackle
*********************************************************************/

#define MU_METAL 16.0

int mean_weight_local_like_grackle(chemistry_data *my_chemistry,
                                   grackle_field_data *my_fields,
                                   gr_float *mu) {

  gr_float number_density = 0.;
  gr_float inv_metal_mol = 1.0 / MU_METAL;
  int i, dim, size = 1;

  for (dim = 0; dim < my_fields->grid_rank; dim++)
    size *= my_fields->grid_dimension[dim];

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  for (i = 0; i < size; i++) {

    if (my_chemistry->primordial_chemistry > 0) {
      number_density =
          0.25 * (my_fields->HeI_density[i] + my_fields->HeII_density[i] +
                  my_fields->HeIII_density[i]) +
          my_fields->HI_density[i] + my_fields->HII_density[i] +
          my_fields->e_density[i];
    }

    /* Add in H2. */

    if (my_chemistry->primordial_chemistry > 1) {
      number_density +=
          my_fields->HM_density[i] +
          0.5 * (my_fields->H2I_density[i] + my_fields->H2II_density[i]);
    }

    if (my_fields->metal_density != NULL) {
      number_density += my_fields->metal_density[i] * inv_metal_mol;
    }

    /* Ignore deuterium. */

    mu[i] = my_fields->density[0] / number_density;
  }

  return SUCCESS;
}

#endif /* MEAN_MOLECULAR_WEIGHT_H */
