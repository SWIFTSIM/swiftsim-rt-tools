#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f 'miniStromgren.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../swiftsim/swift \
    --hydro --threads=4 --stars --external-gravity \
    --feedback --radiation \
    miniStromgren.yml 2>&1 | tee output.log
    # miniStromgrenPetrurbed.yml 2>&1 | tee output.log
    # propagationTest.yml 2>&1 | tee output.log
    # propagationTestPerturbed.yml 2>&1 | tee output.log

# Plot the Stromgren 3D checks.
python3 ./plotSolution.py
