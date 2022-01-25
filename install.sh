#!/bin/bash

#################### DEFINE ENV VARS ####################

# MoBioTools home directory
export MBTHOME=""
export MBTANALYSIS=$MBTHOME/scripts
export PYOVLP=$MBTHOME/src/pyoverlaps
export PATH=$MBTHOME:$PATH
export PATH=$MBTANALYSIS:$PATH
export PATH=$PYOVLP:$PATH


# Molcas
export MOLCAS="" # Add path to molcas build 
export MOLCAS_SCRIPT="" # Add full path of the script 
export MOLCAS_SCRATCH="" 
export OVLPDIR=$PWD

export MOLCAS_MEM=6000
export MOLCAS_MAXMEM=$MOLCAS_MEM
export runmolcas=$MOLCAS_SCRIPT


# Install python modules
python3 setup.py install
