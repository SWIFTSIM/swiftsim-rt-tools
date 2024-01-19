#ifndef SIMULATION_PARAMS_H
#define SIMULATION_PARAMS_H

#include "rt_ngroups.h"
#include "units.h"
#include "parser.h"

struct simulation_params {

  /*! reduced speed of light factor in internal units (units in which simulation will be run) */
  float fc;

  /*! reduced speed of light in internal units (units in which simulation will be run) */
  float c_reduced;

  /*! Particle mass in internal units (units in which simulation will be run) */
  float particle_mass;

  /*! Expected density average in simulation in internal units (units in which simulation will be run) */
  float density_average;

  /*! Expected min for density in internal units (units in which simulation will be run) */
  float density_min;

  /*! Expected max for density in internal units (units in which simulation will be run) */
  float density_max;

  /*! Simulation box size in internal units (units in which simulation will be run) */
  float boxsize;

  /*! Estimate for particle smoothing length in internal units (units in which simulation will be run) */
  float smoothing_length;

  /*! Estimate for minimal radiation energy in sim */
  float rad_energy_min;

  /*! Estimate for maximal radiation energy in sim */
  float rad_energy_max;

  /*! Estimate for average radiation energy in sim */
  float rad_energy_av;



  /*! Are we using constant star emission rates? */
  int use_const_emission_rates;

  /*! Number of particles in the simulation */
  long long npart;



  /*! photon group intervals, in Hz */
  double photon_groups_Hz[RT_NGROUPS];

  /*! Stellar emission rates in internal units (units in which simulation will be run)  */
  double star_emission_rates[RT_NGROUPS];
};


void simulation_params_init(struct simulation_params *p);
void simulation_params_print(struct units *units, struct simulation_params *params);

#endif
