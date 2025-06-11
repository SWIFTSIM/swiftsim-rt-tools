# swiftsim-rt-tools

A collection of Radiative Transfer periphery tools, scripts, and programs for
SWIFT.


## Contents

-   `additionalTests/` Additional test cases which have not made it into the
    main repository, e.g. because they require manual tweaks of the main
    codebase.

-   `benchmarks/` Setups for code benchmarks. Sub-cycling speedup tests,
    scaling tests, etc.

-   `grackle/` - standalone programs using the
    [grackle](https://github.com/grackle-project/grackle-swift) library for specific
    problems.

-   `GEARRTUnitCheck/` - a set of checks for unit compatibility with float
    precision for GEAR RT & Grackle.

-   `HLLRiemannSolverEigenvalues` - compute the Eigenvalues for the HLL Riemann
    solver used (optionally) in the GEAR-RT scheme

-   `integrateCrossSections` - Get average photoionization cross sections,
    averaged over frequency bins.

-   `patches` - patches (`git diff`s) for specific projects/runs/benchmarks.

-   `photonNumberRateToLuminosity` - scripts to convert (ionizing) photon number
    rates into stellar luminosities. Needed e.g. to reproduce standard RT tests.

-   `pyThermochemistry` - some python functions and scripts related to
    thermochemistry. Useful to set up initial conditions, or do some outside
    computations for setups.

-   `readCosmologicalRadiativeTransferComparisonProjectData` - read in data from
    the cosmological radiative transfer comparison project (Iliev et al 2006, 2009)

-   `redshift` - standalone program to let a packet of radiation redshift, either
    accounting for or neglecting an additional redshift interval
