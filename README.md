# swiftsim-rt-tools

A collection of Radiative Transfer periphery tools, scripts, and programs for 
SWIFT.


## Contents

-   `grackle/` - standalone programs using the 
    [grackle](https://github.com/grackle-project/grackle) library for specific 
    problems.

-   `HLLRiemannSolverEigenvalues` - compute the Eigenvalues for the HLL Riemann 
    solver used (optionally) in the GEAR-RT scheme

-   `integrateCrossSections` - Get average photoionization cross sections, 
    averaged over frequency bins.

-   `photonNumberRateToLuminosity` - scripts to convert (ionizing) photon number
    rates into stellar luminosities. Needed e.g. to reproduce standard RT tests.

-   `pyThermochemistry` - some python functions and scripts related to 
    thermochemistry. Useful to set up initial conditions, or do some outside
    computations for setups.


