#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -e IlievTest4ICData.hdf5 ]
then
    echo "Fetching data file for Iliev Test 4 example ..."
    ./getICData.sh
fi

if [ ! -f 'ilievTest4.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../../swiftsim/swift \
    --hydro --threads=4 --stars --external-gravity \
    --feedback --radiation \
    ilievTest4.yml 2>&1 | tee output.log

python3 ./plotSolution.py 9
