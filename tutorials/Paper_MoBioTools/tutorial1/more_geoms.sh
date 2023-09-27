#!/bin/bash
for number in ` seq 200 20 5000` ; do
    sed -i "s/variabless/$number/g" orca.main.in
    main_qminputs.py -i orca.main.in -t orca.tpl.inp
    sed -i "s/$number/variabless/g" orca.main.in
    
done
