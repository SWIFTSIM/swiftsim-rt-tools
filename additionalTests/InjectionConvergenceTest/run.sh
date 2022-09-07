#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail


for res in 16 32 64 128; do
    if [ ! -e glassCube_$res.hdf5 ]
    then
        echo "Fetching initial glass files for Iliev Test 2 example ..."
        ./getGlass.sh
    fi

    if [ ! -f 'ilievTest2-'$res'.hdf5' ]; then
        echo "Generating ICs"
        python3 makeIC.py $res
    fi

    # Run SWIFT with RT
    ../../../swiftsim/swift \
        --hydro --threads=4 --stars --external-gravity \
        --feedback --radiation \
        ilievTest2-$res.yml 2>&1 | tee output-$res.log

done

python3 ./plotSolution.py 10
