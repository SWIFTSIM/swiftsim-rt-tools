#!/bin/bash

set -e

make clean
make
./iliev_test | tee output.log
python3 plot.py
