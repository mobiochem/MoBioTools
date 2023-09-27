#!/bin/bash

currdir=$PWD
n=6

# Generate folder containing xyz files
mkdir xyz_res

for i in {1..5}
do
    # Call main_qminputs.py
    sed "s/closest = 1/closest = "$i"/g" main.inp > main"$i".inp
    main_qminputs.py -i main"$i".inp -t tpl.inp

    # Generate xyz file to visualize geometries
    cd $currdir/geom"$n"
    sed -n -e '/* xyz/,/*/ p' geom"$n".inp > temp; 
    sed -i "1d" temp ; sed -i '$d' temp ; 
    nlines=`wc -l < temp`; 
    sed "1 i $nlines\n" temp > geom"$n".xyz
    cp geom"$n".xyz ../xyz_res/"$i"res-geom"$n".xyz

    cd $currdir
    cd $currdir
    if [ $i != "$n" ]
    then
        rm -r "$i"res-geom"$n"
        mv geom"$n" "$i"res-geom"$n"
    fi


done
