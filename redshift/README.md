A standalone program to test two redshift implementations.

Let a packet of radiation decay with the cosmological expansion of the universe, explicitly accounting for redshifting

## Implementations
In this standalone program, two implementations of explicitly accounting for the redshift integral
$$\int\limits_{\nu_{i,0}}^{\nu_{i,1}}\nu\frac{\partial E_\nu}{\partial\nu}d\nu$$

The first implementation is based on [Finlator et al. (2009)](https://doi.org/10.1111/j.1365-2966.2008.14190.x), and the second is based on the analytic solution of the integral.

## Details
* 10 photon groups spanning $[0, 10\nu_\mathrm{peak}]$
* Initial radiation field assuming blackbody spectrum
* Blackbody spectrum temperature redshifted during time evolution according to $T = T_\mathrm{init} a^{-1}$

This behaviour can be changed in the [include/constants.h](include/constants.h) file. Change the definitions of `RT_NGROUPS`, `SHIFT_BB`, or `FLAT_IC`. 
* `RT_NGROUPS` must be a positive integer
* `SHIFT_BB` must be `1` or `0`
* `FLAT_IC` must be `1` or `0`

## Run
To run the test is straightforward. Three `run.sh` files are shipped. Choose your implementation and run `./run[implementation].sh`. A plot of the energy density decay will be generated and saved under `[implementation].png`. Using the command `./run.sh` runs both implementations and plots both in the same figure under `finlator_analytic_same_plot.png`
