CosmoCoolingTest
-----------------

Run simple cooling test with cosmology enabled.
The output data is intended to provide a reference solution to the
`examples/RadiativeTransferTests/CosmoCoolingTest` example.

A variable `with_cosmo` in `main.c` allows you to turn off cosmological
expansion. In that case, the program will integrate the cooling over
a time frame specified by `a_begin` and `a_end` and your cosmology.
This option is intended to comare results of cosmo and non-cosmo outputs with 
identical ICs run over the identical time frame.

If you compile and run the code with `with_cosmo = 1` and `with_cosmo=2`
separately, you can compare the different results with
`plotComparisonWithNoCosmo.py`

