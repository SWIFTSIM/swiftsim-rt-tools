#!/bin/bash

set -e

make clean
make
./cosmo_cooling_test
./plot.py
