#!/bin/bash


mkdir -p results

for dir in scale_*; do
    rep=${dir#scale_}
    cp -r -u $dir/timesteps.txt results/timesteps-$rep.txt
    cp -r -u $dir/rtsubcycles.txt results/rtsubcycles-$rep.txt
done
