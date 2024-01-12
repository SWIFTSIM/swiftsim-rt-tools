#/bin/bash

# specify number of replications here.
# No reason to do it in a sophisticated way.
replicas="1 2 3 4 5 6 7 8"

#------------------------------------------------------

# fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f swift ]; then
    echo "Didn't find swift executable in this directory"
    exit
fi

# Finally, run SWIFT with RT
workdir=`pwd`
for rep in $replicas; do
    newdir=scale_$rep
    cd $newdir

    jobfile=job-scaling-$rep.sh
    # sbatch $jobfile
    bash $jobfile

    cd $workdir
done
