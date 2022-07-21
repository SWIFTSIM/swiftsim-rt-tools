GEAR-RT: Turning off evolution of hydrodynamics
-----------------------------------------------

This file provides instructions how to modify SWIFT to be able to run GEAR-RT
with Gizmo-MFV hydrodynamics while maintaining a constant density. Several
modifications are required:

1)  Keep particles static. You can either declare all particles in the
    simulation boundary particles at compile time, or you can activate the
    `#GIZMO_FIXED_PARTICLES` macro in `src/const.h`

2)  Disable all particle interactions. E.g. add 

    ```
      dt_therm = 0.f;
      dt_grav = 0.f;
      dt_kick_corr = 0.f;
    ```

    at the top of `src/hydro/Gizmo/hydro.h:hydro_kick_extra()` function.

    You also need to disable the ion mass fluxes in 
    `src/rt/GEAR/rt.h:rt_kick_extra()` by e.g. setting 

    ```
      dt_therm = 0.;
    ```
    at the top of the function (but after the debugging checks...)

3)  While the initial conditions try to generate the exact density requested by 
    the test, the final result may not be exact due to the nature of particle 
    simulations.
    If you want to fix the exact initial density and internal energy, I recommend 
    you do it in the rt_do_thermochemsitry() function.
    At the top, replace the line

    ```
      if (dt == 0.) return;
    ```

    with

    ```
      if (dt == 0.) {
        rt_tchem_set_particle_quantities_for_test(p);
        return;
      }
    ```


    `rt_tchem_set_particle_quantities_for_test()` is a helper function which is
    defined in `src/rt/GEAR/rt_thermochemistry_utils()`.
