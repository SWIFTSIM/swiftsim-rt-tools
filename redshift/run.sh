#!/bin/bash

set -e

make clean
make
./makeIC
./analytic
./finlator

./plot.py -f "ana"
./plot.py -f "fin"
./plotTogether.py
