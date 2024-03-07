#ifndef CHECKS_H
#define CHECKS_H

#include "simulation_params.h"
#include "units.h"

void check_gas_quantities(float density, char *name, float T,
                          const struct simulation_params *params,
                          const struct units *units, int verbose);

void check_grackle_internals(float density, float radiation_energy_density,
                             char *name, float T,
                             const struct simulation_params *params,
                             const struct units *units, int verbose);

void check_radiation_energies(float radEnergy, char *radName, float density,
                              char *name, float T,
                              const struct simulation_params *params,
                              const struct units *units, int verbose);

void check_luminosities(float luminosity, float density, char *name, float T,
                        const struct simulation_params *p,
                        const struct units *units, int verbose);

#endif
