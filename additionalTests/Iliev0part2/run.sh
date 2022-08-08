#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for Iliev Test 0 part 2 example ..."
    ./getGlass.sh
fi

if [ ! -f 'ilievTest0part2.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../../swiftsim/swift \
    --hydro --threads=4 --stars --external-gravity \
    --feedback --radiation \
    ilievTest0part2.yml 2>&1 | tee output.log

# Plot the Stromgren 3D checks.
python3 ./plotSolution.py 101
