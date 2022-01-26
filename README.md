# MoBioTools
A simple yet versatile toolkit to automatically setup quantum mechanics/molecular mechanics (QM/MM) calculations for ensembles of geometries generated from molecular dynamics trajectories.

# Installation Requirements
A python installation with version >= 3.6 is required. Most Linux-based systems have a default python3 installation, otherwise one can use, e.g.: Anaconda (https://www.anaconda.com/products/individual), which includes more than 150 data analysis python libraries, or miniconda (https://docs.conda.io/en/latest/miniconda.html) which is a pocket version of anaconda that includes exclusively python. In either case, new libraries can be installed using the
```
conda
```
command.


Numpy: It comes with the anaconda distribution. Otherwise it can be installed via:

```
conda install numpy
```

```
pip install numpy
```

Cpptraj and pytraj:

If you have an Amber or an AmberTools installation, they should already be present in your $AMBERHOME amber installation directory. Otherwise, cpptraj needs to be installed first. To install pytraj, please visit https://amber-md.github.io/pytraj/latest/installation.html#install. To install cpptraj and pytraj from scratch, visit https://amber-md.github.io/pytraj/latest/installation.html#from-source-code-hard-way-expert-only.

# Installation
Get the MoBioTools source
```
git clone git@github.com:mobiochem/MoBioTools.git
```
then, go to the MoBioTools directory
```
cd MoBioTools
```
set manually the environment variables in the config.sh script, then source it
```
source ./config.sh
```
and finally, carry out the installation itself
```
python3 setup.py install
```


# Usage
The tools consist of two scripts with different functionalities:
1. The input generator (main_qminputs.py, -h option to visualize the arguments it takes), which is compatible with the Gaussian, NWChem and (Open)Molcas quantum mechanics (QM) programs, takes an MD trajectory in amber format and generates an input file for all the snapshots selected by the user. The program takes two input files as arguments: 

- A main input file (-i option) on which the user defines the names of the trajectory file, the topology file, the desired QM program, the QM region to consider among others.

- A template file, on which the specifics of the QM calulation are defined.

Some examples are provided on the examples/inputs folder. 
```
cd examples/inputs
```
For example, in the case of a NWChem calculation, copy the trajectory and the topology files on the examples/inputs folder
```
cp GUA-O2-S5.* nwchem_input/.

cd nwchem_input
```
then execute the main_qminputs.py script as follows:
```
main_qminputs.py -i nwchem.main.in -t nwchem.tpl.inp
```

2. The program to automaticall correct the active space of a set of CASSCF calculations an ensemble of geometries (pyoverlaps.py, -h option to visualize the arguments it takes). 
