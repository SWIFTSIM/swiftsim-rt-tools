#!/bin/bash

set -e

make clean
make
./makeIC
./analytic
./plot.py -f "ana"
