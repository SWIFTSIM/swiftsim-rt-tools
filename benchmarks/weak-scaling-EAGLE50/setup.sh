#!/bin/bash

# specify number of replications here.
# No reason to do it in a sophisticated way.
replicas="1 2 3 4 5 6 7 8"

#------------------------------------------------------

# fail if a subcommand fails
set -e
set -o pipefail


# Create a directory for each subcycle number,
# and create a local copy of the runtime parameters
for REP in $replicas; do

    TOP_LEVEL=$((16*$REP))
    MESH_SIZE=$((512*$REP))
    NMPITASKS=$((4*$REP**3))

    newdir=scale_$REP
    mkdir -p $newdir
    chmod a+w $newdir

    # Copy and adapt parameter files
    ymlfile=./eagle_50_rt_test.yml
    cp $ymlfile $newdir
    # make sure this is the value you also have in the main yml file in this directory
    sed -i "s/REPLICATION/$REP/g" $newdir/$ymlfile
    sed -i "s/MESH_SIZE/$MESH_SIZE/g" $newdir/$ymlfile
    sed -i "s/TOP_LEVEL/$TOP_LEVEL/g" $newdir/$ymlfile

    jobfile=job-scaling-$REP.slm
    cp ./job-EAGLE-weak-scaling.slm $newdir/$jobfile
    sed -i "s/REPLICATION/$REP/g" $newdir/$jobfile
    sed -i "s/NMPITASKS/$NMPITASKS/g" $newdir/$jobfile

done

echo "done setting up"
