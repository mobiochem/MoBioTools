#!/bin/bash
for number in ` seq 199 2 399` ; do
    sed -i "s/variabless/$number/g" gaussian.main.in
    main_qminputs.py -i gaussian.main.in -t gaussian.tpl.inp
    sed -i "s/$number/variabless/g" gaussian.main.in
done
