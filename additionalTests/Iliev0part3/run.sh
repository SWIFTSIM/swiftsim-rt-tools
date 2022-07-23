#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f ./ilievTest0part3.hdf5 ]; then
    echo "creating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../../swiftsim/swift \
    --hydro \
    --threads=4 \
    --verbose=0  \
    --radiation \
    --external-gravity \
    --stars \
    --feedback \
./ilievTest0part3.yml 2>&1 | tee output.log

python3 plotSolution.py
