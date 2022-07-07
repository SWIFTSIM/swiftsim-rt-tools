#!/bin/bash

set -e

make clean
make
./iliev_test
./plot.py
