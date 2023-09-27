#!/bin/bash
apwd=`pwd`

echo '# geom   state   Eexc (eV)  (nm)  f' > ${apwd}/SUMMARY_STEPS.dat

for i in `seq 199 2 399`
do
    for j in `seq 9`
    do
        exc=`grep "Excited State   $j" geom${i}/geom${i}.log | awk '{print $5}'`
        f=`grep "Excited State   $j" geom${i}/geom${i}.log | awk '{print $9}'`
        nm=`grep "Excited State   $j" geom${i}/geom${i}.log | awk '{print $7}'`
        echo "   ${i}   $j    $exc   ${f#*=}   $nm" >> ${apwd}/SUMMARY_STEPS.dat
    done
    exc=`grep "Excited State  10" geom${i}/geom${i}.log | awk '{print $5}'`
    f=`grep "Excited State  10" geom${i}/geom${i}.log | awk '{print $9}'`
    nm=`grep "Excited State  10" geom${i}/geom${i}.log | awk '{print $7}'`
    echo "   ${i}   10    $exc   ${f#*=}   $nm" >> ${apwd}/SUMMARY_STEPS.dat
done
