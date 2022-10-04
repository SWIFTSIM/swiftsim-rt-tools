#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -e glassCube_128.hdf5 ]
then
    echo "Fetching initial glass file for Iliev Test 2 example ..."
    ./getGlass.sh
fi

if [ ! -f 'ilievTest6.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../../swiftsim/swift \
    --hydro --threads=4 --stars --external-gravity \
    --feedback --radiation \
    --steps=1 \
    ilievTest6.yml 2>&1 | tee output.log

python3 ./plotProfiles.py 76
