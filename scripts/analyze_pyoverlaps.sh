#!/bin/bash

# Define variables
# Create results file
if [ -f results.dat ]
then
    rm results.dat
fi

#################### SET MANUALLY RANGE ####################
for i in {0..99}
do
    tail geom"$i"/summary.dat >> results.dat
done

sed -i "/Geom      Last MOs swapped    Corresp. ref. MOs   N iterations        Warnings/d" results.dat
sed -i "1s/^/Geom      Last MOs swapped    Corresp. ref. MOs   N iterations        Warnings\n/" results.dat

# Get conv and nconv
grep "     5     " results.dat > nconv.dat
sed -i "1s/^/Geom      Last MOs swapped    Corresp. ref. MOs   N iterations        Warnings\n/" nconv.dat
sed  "/     5     /d" results.dat > conv.dat

# Get number of geometries recovered and not recovered
num1="$(python3 get_num_geometries.py conv.dat)"
sed -i "1s/^/# $num1 Geometries recovered\n/" conv.dat
num2="$(python3 get_num_geometries.py nconv.dat)"
sed -i "1s/^/# $num2 Geometries not recovered\n/" nconv.dat

