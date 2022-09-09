GEAR-RT: Fix the gas temperature throughout the run
----------------------------------------------------


To keep the gas temperature constant, three things need to be changed.
Firstly, we need to turn off all hydrodynamics. See the corresponding
instruction file on how to do that.

Secondly, we need to turn off the thermochemistry. Look in the 
`rt_do_thermochemistry()` function where the updated gas internal energy
is being set by means of calling `hydro_set_internal_energy(p, u_new);`
and remove that line.

And finally, we need to ensure that the temperature also remains constant 
throughout the run. Even with the thermochemistry and hydrodynamics 
disabled, i.e. when the internal energy remains constant, the temperature
may change because the ionization of the gas changes its mean molecular mass.
So you'll need to add the following lines to `rt_do_thermochemistry()`:


```
  /* Set internal energy to satisfy Temperature requirement for test */
  double mu = rt_tchem_get_mean_molecular_weight(
      p->rt_data.tchem.mass_fraction_HI, p->rt_data.tchem.mass_fraction_HII,
      p->rt_data.tchem.mass_fraction_HeI, p->rt_data.tchem.mass_fraction_HeII,
      p->rt_data.tchem.mass_fraction_HeIII);

  float u_fix = rt_tchem_internal_energy_from_T(
      1.e4, mu, phys_const->const_boltzmann_k, phys_const->const_proton_mass);
  /*  ^---- this is the temperature you set yourself */
  hydro_set_internal_energy(p, u_fix);
```

In the example above, the chosen temperature is 1e4 K.
