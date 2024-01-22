GEAR-RT Units Check
=====================

Run a couple of tests to verify whether your choice of units for a simulation
might be incompatible with the precision limits for the internal thermochemistry
computations.

This program assumes that we're running with Grackle and primordial chemistry =
1, which is the non-equilibrium 6 species network (H0, Hp, He0, Hep, Hepp, e-).

NOTE: these tests don't inlcude cosmology yet.

NOTE: The grackle library is experiencing active development (state 2023). The
API might change in the future. We keep a frozen version forked on
[github](https://github.com/mladenivkovic/grackle-swift). This version is
guaranteed to work with swift (and with this repository).






Instructions
-------------

1)  Generate initial condition data parameters: You need to generate a `.yml`
    file that contains a bunch of initial condition parameters that are necessary
    for this test to be run. See `test/initial_condition_parameter_example.yml` for
    a list of all required parameters. That file is an example parameter yml file
    that the test requires.

    Alternatively, you could run `generate_IC_parameter_file.py` on a SWIFT
    snapshot or IC hdf5 file, which will generate the output for you, provided
    the IC file contains all the required units.

    Note that you can freely edit the resulting IC parameter file
    afterwards, e.g. if you want to test out different min/max particle
    densities you might reach when the simulation evolves.

2)  Define how many photon groups you are planning to run the simulation 
    with at the top of the file `include/rt_ngroups.h`, e.g. 

    ```
    #define RT_NGROUPS 3
    ```

    The default is 3.

3)  Navigate into the `test` directory:
    ```
    $cd test/
    ```

4)  Compile and run the test suite.
    You can compile the code using the provided `Makefile` by running
    ```
    make 
    ```

4.1) Depending on your system setup, you may need to change a
    few variables in the `test/Makefile` file. For example:
    - The directory where grackle is installed (by providing the
      `GRACKLE_SWIFT_ROOT` variable with the correct value)
    - The directory where GSL (gnu scientific library) is installed (by
      providing the `GSL_ROOT` variable with the correct value)
    - The directory where HDF5 (hierarchical data format 5 library) is
      installed (by providing the `HDF5_ROOT` variable with the correct value)
    - The compiler you want to use (by providing the `CC` variable with the
      correct value). The default is gcc.


5)  Run the test suite.
    When compiled successfully in step 4), an executable `test/GEARRT_unit_checks`
    is created. To run the test, you need to run this executable.

    The executable expects 2 command line arguments:
    1)  the .yml parameter file of the ICs or the snapshot you are looking at,
        which you generated in step 1).
    2)  the .yml parameter file you intend to use for your simulation.

    The correct usage would then be something like

    ```
    ./GEARRTUnitCheck path/to/IC_parameters.yml path/to/simulation_parameters.yml
    ```

    Alternatively, you can modify the two arguments hardcoded in the `test/run.sh`
    script. The two `.yml` files have been set up from the `../additionalTests/Iliev5`
    initial conditions, and should work without problems or warnings.


6)  If the tests completed without an issue, there should be 6 `.dat` files written
    as output. Those are results of tests running grackle with the given parameters
    for a cooling and a heating test.
    The results are intended to be visually inspected, as it had occurred in the past
    that weird things kept happening without grackle raising an error.
    To produce plots, simply run the provided scripts
    `tests/plot_cooling_tests.py` and `tests/plot_heating_tests.py`.
    What results you should be expecting is described below.


Results
----------

-   If at any point the code throws an error, it means that most likely your 
    choice of units is not suitable for a run. In any case, it's something
    you definitely need to look into.

-   The code should create two plots: `cooling_test_output.png` and 
    `heating_test_output.png`. Those are results from tests run while using
    grackle where a hot gas parcel is left to cool, and a gas parcel is first
    heated up using radiation and then left to cool, respectively.
    Depending on your initial conditions, the resulting plots may vary 
    strongly. In particular, these tests depends on the number density of the
    particles. The main purpose is to run checks internally while grackle is
    doing its thing, and that the resulting plots should roughly show the 
    characteristics (heating, cooling...) they were set up for. You will most
    likely see strong differences when comparing the results for the average,
    minimal, and maximal density estimate.

-   The test should complete in a matter of seconds. (You can try to establish
    a reference for how fast it should go on your local machine by making a
    direct comparison with the setup for the additional test in 
    swiftsim-rt-tools/additionalTests/Iliev0part3. That finishes in less than
    a second on an intel skylake i7.) If it takes too long, it might indicate 
    that your choice of units is leading to grackle having trouble reaching a 
    convergent solution. Try changing your units and see whether that can make 
    the test finish in seconds.


