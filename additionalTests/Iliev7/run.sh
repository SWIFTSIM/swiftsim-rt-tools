#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f 'ilievTest7.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../../swiftsim/swift \
    --hydro --threads=4 --stars --external-gravity \
    --feedback --radiation \
    --steps=1 \
    ilievTest7.yml 2>&1 | tee output.log
    # ilievTest7UnequalMasses.yml 2>&1 | tee output.log

python3 ./plotSlices.py 51
