#!/bin/bash


mkdir -p results

for dir in stromgren_test_*; do
    cp -r -u $dir/timing* results
    nsub=${dir#stromgren_test_}
    cp -r -u $dir/timesteps.txt results/timesteps-$nsub.txt
done
