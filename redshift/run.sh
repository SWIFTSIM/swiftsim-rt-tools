#!/bin/bash

set -e

make clean
make
./cosmo_redshift_test
./plot.py
