#ifndef PHOTON_INTERACTION_RATES_H
#define PHOTON_INTERACTION_RATES_H

#include "constants.h"
#include "cross_sections.h"

/* ------------------------------------------------------------- */
/* Get the heating and ionization rates from photon interations  */
/* ------------------------------------------------------------- */

/**
 * @brief compute the heating, ionization, and dissassociation rates
 * for the particle radiation field as needed by grackle, and the
 * net absorption/emission rates for each photon group
 *
 * NOTE: this assumes that all quantities are in cgs
 *
 * @param radiation_energy_density energy densities of each frequency group
 * @param species_densities the physical densities of all traced species
 * @param energy_weighted_cross_sections: energy weighted averaged cross
 *sections.
 * @param number_weighted_cross_sections: number weighted averaged cross
 *sections.
 * @param mean_energy: mean photon energy in each photon frequency bin.
 * @param time_units: internal time units conversion factor to cgs.
 *                    val * time_units = val in cgs
 * @param rates (return) Interaction rates for grackle. [0]: heating rate.
 * [1]: HI ionization. [2]: HeI ionization. [3]: HeII ionization.
 * [4]: H2 dissociation.
 **/
void get_interaction_rates(double radiation_energy_density[RT_NGROUPS],
                           gr_float species_densities[6],
                           double **energy_weighted_cross_sections,
                           double **number_weighted_cross_sections,
                           double mean_energy[RT_NGROUPS], double time_units,
                           gr_float rates[5]) {

  rates[0] = 0.; /* Needs to be in [erg / s / cm^3 / nHI] for grackle. */
  rates[1] = 0.; /* [1 / time_units] */
  rates[2] = 0.; /* [1 / time_units] */
  rates[3] = 0.; /* [1 / time_units] */
  rates[4] = 0.; /* [1 / time_units] */

  /* "copy" ionization energies from cross section parameters */
  struct photoion_cs_parameters cs_params_cgs = init_photoion_cs_params_cgs();
  const double *E_ion = cs_params_cgs.E_ion;

  /* First, get species number densities and number densities */
  /* in units of neutral hydrogen number density. */
  double species_number_densities[RT_NIONIZING_SPECIES]; /* in cm^-3 */
  /* neutral hydrogen */
  species_number_densities[0] = species_densities[0] / const_mh;
  /* neutral helium */
  species_number_densities[1] = 0.25 * species_densities[2] / const_mh;
  /* singly ionized helium */
  species_number_densities[2] = 0.25 * species_densities[3] / const_mh;

  /* store photoionization rate for each species here */
  double ionization_rates_by_species[RT_NIONIZING_SPECIES];
  for (int spec = 0; spec < RT_NIONIZING_SPECIES; spec++) {
    ionization_rates_by_species[spec] = 0.;
  }

  for (int group = 0; group < RT_NGROUPS; group++) {

    /* Sum results for this group over all species */
    /* double heating_rate_group_nHI = 0.; */
    double heating_rate_group = 0.;
    const double Emean = mean_energy[group];
    const double Eic = radiation_energy_density[group] * const_speed_light_c;
    double Nic = 0.;
    if (Emean > 0.)
      Nic = Eic / Emean;

    for (int spec = 0; spec < RT_NIONIZING_SPECIES; spec++) {
      /* Note: the cross sections are in cgs. */
      const double cse = energy_weighted_cross_sections[group][spec];
      const double csn = number_weighted_cross_sections[group][spec];
      const double n_spec = species_number_densities[spec];

      heating_rate_group += (cse * Emean - E_ion[spec] * csn) * n_spec;
      ionization_rates_by_species[spec] += csn * Nic;
    }

    /* rates[0] += heating_rate_group * Eic / nHI; */
    rates[0] += heating_rate_group * Nic;
  }

  /* Unit conversions for grackle */
  for (int spec = 0; spec < RT_NIONIZING_SPECIES; spec++) {
    /* Grackle wants them in 1/internal_time_units */
    ionization_rates_by_species[spec] /= (1. / time_units);
  }
  /* Grackle wants heating rate in units of / nHI_cgs */
  rates[0] /= species_number_densities[0];

  /* We're done. Write the results in correct place */
  rates[1] = ionization_rates_by_species[0];
  rates[2] = ionization_rates_by_species[1];
  rates[3] = ionization_rates_by_species[2];
  /* rates[4] = skipped for now */
}

#endif
