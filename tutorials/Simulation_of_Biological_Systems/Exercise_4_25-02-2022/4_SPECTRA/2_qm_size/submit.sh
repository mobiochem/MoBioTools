#!/bin/bash
n=6

for i in {1..5}
do
    cd "$i"res-geom"$n"
    orca geom"$n".inp > geom"$n".inp.out
    cd ..
done
