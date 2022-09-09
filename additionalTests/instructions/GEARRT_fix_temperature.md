GEAR-RT: Fix the gas temperature throughout the run
----------------------------------------------------


To keep the gas temperature constant, two things need to be changed.
Firstly, we need to turn off all hydrodynamics. See the corresponding
instruction file on how to do that.

Secondly, we need to turn off the thermochemistry. Look in the 
`rt_do_thermochemistry()` function where the updated gas internal energy
is being set by means of calling `hydro_set_internal_energy(p, u_new);`
and remove that line.
