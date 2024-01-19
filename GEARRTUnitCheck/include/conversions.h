#ifndef CONVERSIONS_H
#define CONVERSIONS_H

#include "simulation_params.h"
#include "units.h"

float conversions_estimate_dt(const struct simulation_params *p);
float conversions_radiation_energy_density_from_luminosity(double luminosity, const struct simulation_params *p, const struct units* units);


#endif
