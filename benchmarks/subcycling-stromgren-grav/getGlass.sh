#!/bin/bash

# read in resolution
read line < resolution;
res=$line;

wget http://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/glassCube_"$res".hdf5

