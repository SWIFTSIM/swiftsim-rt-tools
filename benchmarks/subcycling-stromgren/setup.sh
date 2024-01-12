#!/bin/bash

# specify number of subcycles here. 
# No reason to do it in a sophisticated way.
nsubcycles="1 2 4 8 16 32 64 128 256"

#------------------------------------------------------

# fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f resolution ]; then
    echo "Didn't find 'resolution' file"
    exit
fi

# read in resolution
read line < resolution;
res=$line;

# which glass file to fetch?
glassfile="glassCube_"$res".hdf5"

# fetch it now
if [ ! -f "$glassfile" ]
then
    echo "Fetching initial glass file for Strömgen Sphere 3D example ..."
    ./getGlass.sh
fi

# generate ICs?
# Uncomment this if you want the hydrogen-only ICs too
# if [ ! -f "stromgrenSphere-3D-$res.hdf5" ]; then
#     echo "Generating ICs"
#     python3 deprecated/makeIC.py
#     # srun -n 1 --partition p4 python3 makeIC.py
# else
#     echo found H IC file
# fi
if [ ! -f "stromgrenSphere-3D-HHe-$res.hdf5" ]; then
    echo "Generating ICs"
    python3 makeIC_HHe.py
    # srun -n 1 --partition p4 python3 makeIC_HHe.py
else
    echo found He IC file
fi

# Create a directory for each subcycle number,
# and create a local copy of the runtime parameters
for n in $nsubcycles; do
    newdir=stromgren_test_$n
    mkdir -p $newdir

    # for run in MF MFHHe singlebin; do
    for run in MFHHe; do
        # Copy and adapt parameter files
        ymlfile=./stromgrenSphere-3D-$run-$res.yml 
        cp $ymlfile $newdir
        # make sure this is the value you also have in the main yml file in this directory 
        sed -i "s/max_nr_rt_subcycles: 1/max_nr_rt_subcycles: $n/" $newdir/$ymlfile
        sed -i "s;file_name: ../stromgrenSphere-3D.hdf5;file_name: ../stromgrenSphere-3D-$res.hdf5;" $newdir/$ymlfile
        sed -i "s;file_name: ../stromgrenSphere-3D-HHe.hdf5;file_name: ../stromgrenSphere-3D-HHe-$res.hdf5;" $newdir/$ymlfile

        # without gravity
        jobfile=job-stromgren-$run-$res.slm  
        cp ./job-stromgren-subcycling.slm $newdir/$jobfile
        sed -i "s/#SBATCH --job-name=StromSS/#SBATCH --job-name=Sub$run-$n-$res/" $newdir/$jobfile
        sed -i "s;paramfile.yml;$ymlfile;" $newdir/$jobfile
        sed -i "s;timingfile;timing-$res-$run-$n;" $newdir/$jobfile

        # including gravity
        # jobfile=job-stromgren-$run-$res-grav.slm
        # cp ./job-stromgren-subcycling.slm $newdir/$jobfile
        # sed -i "s/#SBATCH --job-name=StromSS/#SBATCH --job-name=Sub$run-$n-$res/" $newdir/$jobfile
        # sed -i "s/--external-gravity/--self-gravity/" $newdir/$jobfile
        # sed -i "s;paramfile.yml;$ymlfile;" $newdir/$jobfile
        # sed -i "s;timingfile;timing-$res-$run-$n-grav;" $newdir/$jobfile
    done

    # jobfile=job-stromgren-singlebin-$res.slm
    # sed -i "s;../swift ;../swift_singlebin ;" $newdir/$jobfile
done

echo "done setting up"
