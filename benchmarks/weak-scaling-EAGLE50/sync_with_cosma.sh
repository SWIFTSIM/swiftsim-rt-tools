#!/bin/bash

#------------------------------------------------
# rsync my scripts to/from given destination
# using ssh.
#------------------------------------------------




function errormsg(){
    echo "I need you to tell me in which direction to sync."
    echo "Usage: "
    echo "   sync_with_cosma.sh local"
    echo "        to sync TO local machine (e.g. FROM cosma machine)"
    echo "   sync_with_cosma.sh cosma"
    echo "        to sync TO cosma machine"
    exit
}



if [[ `hostname` == *"cosma"* ]] || [[ `hostname` == *"dataweb"* ]]; then
    echo DO NOT USE THIS SCRIPT FROM COSMA OR FROM DATAWEB SERVER.
    echo USE IT FROM YOUR LOCAL MACHINE.
    exit
fi




if [[ "$#" < 1 ]]; then
    errormsg;
else
    case $1 in

        cosma)
            echo "Sending to COSMA"
            UNAME=dc-ivko1
            HOST=dataweb.cosma.dur.ac.uk
            SRC=$HOME/coding/simulations/ready-to-go/cosma8/swift-rt-weak-scaling-EAGLE
            DESTDIR=/cosma/home/do009/dc-ivko1/runs_swift/rt-weak-scaling-EAGLE50
            DEST="$UNAME"@"$HOST":"$DESTDIR"
        ;;

        local)
            echo "syncing to local machine"
            UNAME=dc-ivko1
            HOST=dataweb.cosma.dur.ac.uk
            SRCDIR=/cosma/home/do009/dc-ivko1/runs_swift/rt-weak-scaling-EAGLE50
            SRC="$UNAME"@"$HOST":"$SRCDIR"
            DEST=$HOME/coding/simulations/ready-to-go/cosma8/swift-rt-weak-scaling-EAGLE
        ;;

        *)
            errormsg;
        ;;

    esac

fi



rsync                   \
    --verbose           \
    --recursive         \
    --update            \
    --links             \
    --hard-links        \
    --executability     \
    --perms             \
    --times             \
    --human-readable    \
    -e ssh              \
    --exclude .git      \
    --exclude '*hdf5'   \
    --exclude 'eagle_test_*'  \
    "$SRC" \
    "$DEST"
