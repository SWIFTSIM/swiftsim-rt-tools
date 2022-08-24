GEAR-RT: Get ionization front radius each step
----------------------------------------------------

For some idealized tests, the current position and velocity of the I-front are
required outputs. For higher resolution runs, it could be too costly to write an
output every handful of steps, while interpolating the results from only a handful
of snapshots is not accurate enough (in particular for the velocities.)


Solution: Collect and write out the current position of the I-front each step.

Note: The solution documented below only works for runs without MPI, and assumes
that there is only 1 star.


Instructions
-----------------

1)  Add a file pointer to the `struct rt_properties`, e.g.

    ```
    struct rt_props {

      /* lots of stuff */

      FILE* r_ifront_fp;
    };

    ```


2) Open the file in `rt_props_init()` and write a header, e.g.:

    ```

    __attribute__((always_inline)) INLINE static void rt_props_init(..) {
    
      /* lots of stuff */

      rtp->r_ifront_fp = fopen("ionization_front_radius.txt", "w");
      fprintf(rtp->r_ifront_fp, "# length units: %.12e\n",
              units_cgs_conversion_factor(us, UNIT_CONV_LENGTH));
      fprintf(rtp->r_ifront_fp, "# time units: %.12e\n",
              units_cgs_conversion_factor(us, UNIT_CONV_TIME));
      fprintf(rtp->r_ifront_fp, "# %14s %14s", "time [IU]", "radius [IU]\n");

    }
    ```

3)  Add the following code to the top of the file `engine.c` (the functions must
    be defined before they are referenced):

```
    struct histogram_data {
      int nbins;
      int *count;
      double *hist;
      double sx[3];
      double boxsize;
    };

    void rt_histogram_xHI_mapper(void *restrict map_data, int count,
                                 void *restrict extra_data) {

      struct part *restrict parts = (struct part *)map_data;
      struct histogram_data *restrict hist_data =
          (struct histogram_data *)extra_data;

      int *count_loc = malloc(hist_data->nbins * sizeof(int));
      double *hist_loc = malloc(hist_data->nbins * sizeof(double));

      for (int i = 0; i < hist_data->nbins; i++) {
        count_loc[i] = 0;
        hist_loc[i] = 0.;
      }

      double dr = 0.5 * hist_data->boxsize / hist_data->nbins;

      for (int k = 0; k < count; k++) {

        struct part *restrict p = &parts[k];
        double xHI = (double)p->rt_data.tchem.mass_fraction_HI /
                     ((double)p->rt_data.tchem.mass_fraction_HI +
                      (double)p->rt_data.tchem.mass_fraction_HII);

        /* Ignore periodicity */
        double dx = p->x[0] - hist_data->sx[0];
        double dy = p->x[1] - hist_data->sx[1];
        double dz = p->x[2] - hist_data->sx[2];
        double d = sqrt(dx * dx + dy * dy + dz * dz);
        if (d > 0.5 * hist_data->boxsize) continue;

        int index = floor(d / dr);
        hist_loc[index] += xHI;
        count_loc[index] += 1;
      }

      for (int i = 0; i < hist_data->nbins; i++) {
        atomic_add_d(&(hist_data->hist[i]), hist_loc[i]);
        atomic_add(&(hist_data->count[i]), count_loc[i]);
      }
    }

    void rt_get_Ifront_radius(struct engine *e) {

      struct rt_props *rt_props = e->rt_props;

      /* Get current time, valid for both main steps and subcycles */
      double time = e->ti_current_subcycle * e->time_base + e->time_begin;

      const int nbins = 10;
      int *count = malloc(nbins * sizeof(int));
      double *hist = malloc(nbins * sizeof(double));
      double *profile = malloc(nbins * sizeof(double));

      for (int i = 0; i < nbins; i++) {
        count[i] = 0;
        hist[i] = 0.;
        profile[i] = 0.;
      }

      struct histogram_data hist_data;
      hist_data.count = count;
      hist_data.hist = hist;
      hist_data.nbins = nbins;
      hist_data.boxsize = e->s->dim[0];

      /* Assume only 1 star in simulation */
      struct spart s = e->s->sparts[0];
      hist_data.sx[0] = s.x[0];
      hist_data.sx[1] = s.x[1];
      hist_data.sx[2] = s.x[2];

      /* threadpool */
      if (e->s->nr_parts > 0)
        threadpool_map(&e->threadpool, rt_histogram_xHI_mapper, e->s->parts,
                       e->s->nr_parts, sizeof(struct part),
                       threadpool_auto_chunk_size, /*extra_data=*/&hist_data);

      /* Compute the actual profile */
      for (int i = 0; i < nbins; i++) {
        double res = 0.;
        if (count[i] > 0) res = hist[i] / (double)count[i];
        profile[i] = res;
      }

      /* Now find the ionization front radius */
      double rI = -1.;
      double dr = (0.5 * e->s->dim[0]) / (double)nbins;

      for (int i = 0; i < nbins; i++) {

        /* I-front is defined as 50% neutral fraction */
        if (profile[i] >= 0.5) {
          if (i == 0) {
            rI = 0.;
          } else {
            /* Interpolate exact position */
            double dx = profile[i] - profile[i - 1];
            double a = dx / dr;
            /* assume bin values are represented at center.
             * this bin, where xH > 0.5, has center at (i + 0.5) * dr
             * then left interpl. point is at profile[i-1], (i - 0.5) * dr */
            double r1 = ((double)i - 0.5) * dr;
            double b = profile[i - 1] - a * r1;
            rI = (0.5 - b) / a;
          }
          break;
        }
      }

      fprintf(rt_props->r_ifront_fp, "%12.6e %12.6e\n", time, rI);
    }
```

4) Call the actual functions *after* each step. E.g. in `engine_step()`:

    ```

    int engine_step(struct engine *e) {

      /* lots of stuff */

      rt_get_Ifront_radius(e);

      /* Collect information about the next time-step */
      engine_collect_end_of_step(e, 1);

      /* lots of stuff */
      
    }
    ```

    If you are running with subcycling, you'll probably want to write to the
    file after each sub-cycle as well. So you'll need to add the same function
    call to `engine_run_rt_sub_cycles()`:


    ```

    void engine_run_rt_sub_cycles(struct engine *e) {

      /* Lots of stuff */

      for (int sub_cycle = 1; sub_cycle < nr_rt_cycles; ++sub_cycle) {

        /* Lots of stuff */

        rt_get_Ifront_radius(e);

        /* Collect number of updates and print */
        engine_collect_end_of_sub_cycle(e);

        /* Lots of stuff */
      }
    }

    ```
