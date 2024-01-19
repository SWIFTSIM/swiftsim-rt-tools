#include "simulation_params.h"
#include "constants.h"
#include "conversions.h"
#include "error.h"

/**
 * Initialize an empty simulation params struct.
 **/
void simulation_params_init(struct simulation_params *p) {

  p->fc = 0.;
  p->c_reduced = 0.;

  p->particle_mass = 0.f;
  p->density_average = 0.f;
  p->density_min = 0.f;
  p->density_max = 0.f;
  p->boxsize = 0.f;
  p->smoothing_length = 0.f;

  p->rad_energy_min = 0.f;
  p->rad_energy_max = 0.f;
  p->rad_energy_av = 0.f;

  p->use_const_emission_rates = 0;
  p->npart = 0ll;

  for (int g = 0; g < RT_NGROUPS; g++) {
    p->photon_groups_Hz[g] = 0.;
    p->star_emission_rates[g] = 0.;
  }
}

/**
 * @brief print out the used parameters for a visual inspection
 **/
void simulation_params_print(struct units *units,
                             struct simulation_params *params) {

  message("Units: [in cgs]");
  message("%22s: %.6e", "mass units", units->mass_units);
  message("%22s: %.6e", "length units", units->length_units);
  message("%22s: %.6e", "time units", units->time_units);
  message("%22s: %.6e", "density units", units->density_units);
  message("%22s: %.6e", "velocity units", units->velocity_units);
  message("%22s: %.6e", "temperature units", units->temperature_units);
  message("%22s: %.6e", "energy units", units->energy_units);
  message("%22s: %.6e", "internal energy units", units->internal_energy_units);

  message("");
  message("Variables [internal units]");
  message("%22s: %.6e", "particle mass", params->particle_mass);
  message("%22s: %.6e", "density_average", params->density_average);
  message("%22s: %.6e", "density_min", params->density_min);
  message("%22s: %.6e", "density_max", params->density_max);
  message("%22s: %.6e", "boxsize", params->boxsize);
  message("%22s: %.6e", "smoothing length", params->smoothing_length);
  message("%22s: %.6e", "approx dt [internal units]",
          conversions_estimate_dt(params));
  message("%22s: %.6e", "approx dt [s]             ",
          conversions_estimate_dt(params) / units->time_units);
  message("%22s: %.6e", "approx dt [kyr]           ",
          conversions_estimate_dt(params) * units->time_units / const_yr *
              1e-3);
}
