import numpy as np

# Define conversion factor between kcal/mol and eV
conv = 1/23.0605

# Define variables containing input files
ox = "oxidized/run/energies.dat"
red = "reduced/run/energies.dat"

# read the energies for the oxidized and for the reduced species
# each set of energies will be contained in one array

Eox = np.loadtxt(ox, usecols = (1))
Ered = np.loadtxt(red, usecols = (1))

# Print either of the arrays

# Compute averages
av_ox = np.average(Eox)
av_red = np.average(Ered)

# Compute Gibbs Free Energy of reduction
# Also multiply by conv
G = conv * (av_red - av_ox)
print("Gibbs free energy of reduction = {} eV".format(G))
