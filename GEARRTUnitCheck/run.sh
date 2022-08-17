#!/bin/bash

# quit if something fails
set -e

make clean
make
./GEARRT_unit_checks
python3 plot_cooling_tests.py
python3 plot_heating_tests.py
