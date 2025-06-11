#!/bin/bash

set -e

make clean
make
./makeIC
./finlator
./plot.py -f "fin"
