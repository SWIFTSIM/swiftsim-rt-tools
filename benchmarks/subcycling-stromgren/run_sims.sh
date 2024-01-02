#/bin/bash

# specify number of subcycles here. 
# No reason to do it in a sophisticated way.
# nsubcycles="1 2 4 8 16 32"
nsubcycles="1 2 4 8 16 32 64 128"

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

if [ ! -f swift ]; then
    echo "Didn't find swift executable in this directory"
    exit
fi

# if [ ! -f swift_singlebin ]; then
#     echo "Didn't find swift executable in this directory"
#     exit
# fi

# Finally, run SWIFT with RT
workdir=`pwd`
for n in $nsubcycles; do
    newdir=stromgren_test_$n
    cd $newdir

    # for run in MF MFHHe singlebin; do
    # for run in MF MFHHe; do
    for run in MFHHe; do
        jobfile=job-stromgren-$run-$res.slm  
        sbatch $jobfile
        jobfile=job-stromgren-$run-$res-grav.slm  
        sbatch $jobfile
    done

    cd $workdir
done
