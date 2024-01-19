#ifndef UNITS_H
#define UNITS_H

/*! A struct to store units. */
struct units {
  double mass_units;
  double time_units;
  double length_units;
  double density_units;
  double velocity_units;
  double temperature_units;
  double internal_energy_units;
  double energy_units;
  double energy_density_units;
  double power_units;
};

void units_init(struct units *units);
void units_get_internal_units(struct units *u);

#endif
