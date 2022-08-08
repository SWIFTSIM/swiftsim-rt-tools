GEAR-RT: (Temporarily) Keep Radiation Fields Fixed
----------------------------------------------------

These instructions describe how to maintain a constant radiation field for 
GEAR-RT and then turn it off after the simulation has reached some time.

To fix the radiation fields depending on the current simulation time, add a 
`const double time` parameter to `rt_do_thermochemistry()`. This function is 
particularly good because it already has all the units etc passed as parameters, 
and is the last step right before we use the radiation fields. You get the time 
in internal units via `e->time`, or in case of sub-cycling (discouraged) 
`e->ti_current_subcycle * e->time_base`, which is available when the 
thermochemistry task is called.  Then add the following line to 
`rt_do_thermochemistry()` at the top (important: BEFORE the call to 
`rt_tchem_get_interaction_rates()`)

```
  rt_tchem_set_particle_radiation_field(p, time, us);
```

Secondly, if you are using the cooling time as an upper limit for the particle 
time step, then you will want the same radiation fields to be present when the
cooling time is being computed. To achieve that, either remove the lines in
`rt_do_thermochemistry()` where the radiation fields are updated, or
alternatively just add the same line as above at the bottom of the function.

The function `rt_tchem_set_particle_radiation_field()` is a helper function which is
defined in `src/rt/GEAR/rt_thermochemistry_utils()`. You can set the radiation
field in there, and have access to the time so you may turn it on or off
whenever you like.
