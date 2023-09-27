#!/bin/bash

# Minimization
sander -O -i min.in -o min.out -c MOL-red.inpcrd -p MOL-red.prmtop -r min.rst7 

# Heating
sander -O -i heat.in -o heat.out -c min.rst7 -p MOL-red.prmtop -r heat.rst7 -x heat.mdcrd

# Classical MD
sander -O -i md.in -o md.out -c heat.rst7 -p MOL-red.prmtop -r md.rst7 -x md.mdcrd

# QM/MM MD
sander -O -i qmmm.in -o qmmm.out -c md.rst7 -p MOL-red.prmtop -r qmmm.rst7 -x qmmm.mdcrd
#sander -O -i qmmm.in -o qmmm.out -c MOL.inpcrd -p MOL.prmtop -r qmmm.rst7 -x qmmm.mdcrd
