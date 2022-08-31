#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -e glassCube_128.hdf5 ]
then
    echo "Fetching initial glass file for Iliev Test 2 example ..."
    ./getGlass.sh
fi

if [ ! -f 'ilievTest2.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../../swiftsim/swift \
    --hydro --threads=4 --stars --external-gravity \
    --feedback --radiation \
    ilievTest2.yml 2>&1 | tee output.log

# Plot the Stromgren 3D checks.
python3 ./plotSolution.py 10
