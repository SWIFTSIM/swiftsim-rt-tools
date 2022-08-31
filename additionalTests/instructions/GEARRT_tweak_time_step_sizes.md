GEAR-RT: Manually Tweak Time Step Sizes
----------------------------------------------------

For tests like the Iliev 0 Part 3 test, the reference solutions require the test
to start with very small time steps and end at comparatively very large end
times, leading to millions or billions of required time steps. Simultaneously,
the time steps are allowed to increase over the simulation time. So this file
documents how to make the RT time steps increase over time.

This hasn't been tested with subcycling turned on. I don't recommend you use it.


First, add a time step size to the RT particle struct:

in `rt_struct.h`:

```
    /* Additional RT data in hydro particle struct */
    struct rt_part_data {
        
      /* lots of stuff */

      double dt;

      /* lots of stuff */
    }
```


Then initialize it to some value you want in `rt_first_init_part()` in `rt.h`:


```
    void rt_first_init_part(
        struct part* restrict p, const struct cosmology* cosmo,
        const struct rt_props* restrict rt_props) {

      /* stuff */

      /* Obviously the initial value needs to be in internal units already. */
      p->rt_data.dt = 1.597603e-10;
    }
```


Change what `rt_compute_timestep()` returns in `rt.h`:


```
    static float rt_compute_timestep(
        const struct part* restrict p, const struct xpart* restrict xp,
        struct rt_props* rt_props, const struct cosmology* restrict cosmo,
        const struct hydro_props* hydro_props,
        const struct phys_const* restrict phys_const,
        const struct unit_system* restrict us) {

      /* stuff */

      /* replace current return value with this line */
      return min3(p->rt_data.dt, dt, fabsf(dt_cool));
    }
```


Finally, increase the value each step. I recommend the modification to be done
in `rt_do_thermochemistry()`, e.g. like this:


```
    static void rt_do_thermochemistry(
        struct part* restrict p, struct xpart* restrict xp,
        struct rt_props* rt_props, const struct cosmology* restrict cosmo,
        const struct hydro_props* hydro_props,
        const struct phys_const* restrict phys_const,
        const struct unit_system* restrict us, const double dt, const double time) {

      /* Lots of stuff */

      p->rt_data.dt *= 1.01;

      /* Lots of stuff */

    }
```
