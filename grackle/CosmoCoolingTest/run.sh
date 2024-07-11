#!/bin/bash

set -e

make clean
make
./cosmo_cooling_test | tee output_cosmo.log
./cosmo_cooling_test_physical | tee output_physical.log
./plotComparisonWithPhysical.py
