#!/bin/bash

for i in {0..9}
do
    cd geom"$i"
    orca geom"$i".inp > geom"$i".inp.out
    cd ..
done
