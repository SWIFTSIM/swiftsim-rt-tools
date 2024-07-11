CosmoCoolingTest
-----------------

Run simple cooling test with cosmology enabled.
The output data is intended to provide a reference solution to the
`examples/RadiativeTransferTests/CosmoCoolingTest` example.

A variable `with_cosmo` in `main.c` allows you to turn off cosmological
expansion. In that case, the program will integrate the cooling over
a time frame specified by `a_begin` and `a_end` and your cosmology.
This option is intended to comare results of cosmo and non-cosmo outputs with 
identical ICs run over the identical time frame to identify initial conditions
with visible differences for tests.

If you compile and run the code with `with_cosmo = 1` and `with_cosmo=0`
separately, you can compare the different results with
`plotComparisonWithNoCosmo.py`

A second program has been added, stored in `main_physical.c`. It
does the same as `main.c`, in particular it accounts for the expansion of
the universe, but grackle is given all quantities in proper units and
doesn't use grackle's cosmological units.

Compare results obtained with `main.c` and `main_physical.c` using
`plotComparisonWithPhysical.py`
