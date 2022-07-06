#!/bin/bash

set -e

make clean
make
./cooling_test
./plot.py
