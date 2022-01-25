import numpy as np
from math import pi, exp
from time import process_time
import sys
import pickle


def save_matrix(M, prefix, suffix, symmetric = True):
    with open(prefix + suffix, "w") as f:
        fmt = "{:>7d} {:> 0.6e} {:> 0.6e} {:> 0.6e} {:> 0.6e} {:> 0.6e}"
        for i in range(int(len(M[0])/5)): # Loop over columns
            ind = np.array([k + 5*i for k in range(0,5)]) 
            f.write("{:>17d}{:14d}{:14d}{:14d}{:14d}".format(*(ind + 1)) + "\n")
            if(symmetric):
                inf_inner = 5*i
            else:
                inf_inner = 0
            for j in range(inf_inner, len(M)): # Loop over rows
                f.write(str(fmt.format(j+1, M[j][ind[0]], M[j][ind[1]], 
                    M[j][ind[2]], M[j][ind[3]], M[j][ind[4]])).replace("e", "D") + "\n")

