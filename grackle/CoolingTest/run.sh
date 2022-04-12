#!/bin/bash

set -e

make clean
make
./rt
./plot.py
