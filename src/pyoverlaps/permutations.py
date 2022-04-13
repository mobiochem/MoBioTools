#!/usr/bin/env python3
import numpy as np
from time import process_time

# Non-standard (i.e. non-pyscf library) Basis sets
ATOMNUM = {1:"H", 3: "Li", 4: "Be", 6: "C", 7: "N", 8: "O", 9: "F",\
           11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S",\
           17: "Cl", 19: "K", 20: "Ca", 21: "Sc", 22:"Ti", 23:"V",\
           24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu",\
           30: "Zn", 35: "Br", 46: "Pd", 47: "Ag", 78: "Pt", 79: "Au"}

def readcrd(infile):
    """Read atom coordinates from any file containing
    the coordinate format: atom x y z.
    Returns a list of lists in the internal atom format of
    pyscf: [[atom1, (x1, y1, z1)], [atom2, (x2, y2, z2)], ...]
    """
    coords = []

    with open(infile, "r") as f:
        lines = f.read().splitlines()

    for cnt, iline in enumerate(lines):
        row = iline.split()
        if((len(row) >= 4) and (row[0] in ATOMNUM.keys()\
                and (row[1] not in ATOMNUM.keys()))):
            atom = row[0]
            xyz = tuple(row[-3:])
            coords.append([atom, xyz])
    return(coords)
        


def perm(arg, ATOMNUM):
    """Permute elements on a mol._bas section corresponding
    to an atom entry if they correspond to a SP shell"""
    arglist = np.array(arg)
    if(ATOMNUM<=10):
#        permid = arglist[:2] + arglist[[3,2]] + arglist[4:]
        permid = np.concatenate((arglist[:2], arglist[[3,2]], arglist[4:]))
    elif(ATOMNUM<=18):
        print(arglist[[4,2,5,3,6]])
#        permid = arglist[:2] + arglist[[4,2,5,3,6]] + arglist[7]
        permid = np.concatenate((arglist[:2], arglist[[4,2,5,3,6]], arglist[7:]))
    return(list(permid))

def get_d_indices(mol):
    """mol = Mole() object withh _bas attribute. Returns list
    with d orbital indices"""
    count = 0
    idlist = []
    for cnt, iarr in enumerate(mol._bas):
        amom = iarr[1]
        if(amom == 0):
            count += 1
        elif(amom == 1):
            count += 3
        elif(amom == 2):
            idlist += list(range(count, count + 5))
            count +=5
        elif(amom ==3):
            count +=7
    return(idlist)

def perm_d_indices(mol):
    """Similar to get_d_indices, but it also performs the index swapping
    so that it corresponds to the gaussian09 ordering. It returns the
    full ipermuted index list.

    UPDATE: Extend loop to f functions, do not swap them
    """
    count = 0
    idlist = []
    for cnt, iarr in enumerate(mol._bas):
        amom = iarr[1]
        if(amom == 0):
            idlist += [count]
            count += 1
        elif(amom == 1):
            idlist += list(range(count, count + 3))
            count += 3
        elif(amom == 2):
            temp = np.arange(count, count + 5)
            idlist += list(temp[[2, 3, 1, 4, 0]])
#            idlist += list(temp[[4, 2, 0, 1, 3]])
            count +=5
        elif(amom==3):
            temp = np.arange(count, count + 7)
            idlist += list(temp[[3, 4, 2, 5, 1, 6, 0]])
#            idlist += list(temp[[6, 4, 2, 0, 1, 3, 5]])
#            idlist += list(range(count, count + 7))
            count +=7

    return(idlist)

def perm_matrix(M, indices):
    """Permute matrix using indexing on the indices argument"""
    return(np.transpose(M[indices])[indices])

def perm_mixed_matrix(M, ind1, ind2):
    """Permute a non-symmetric matrix using indexing on the indices argument"""
    dim = len(M)
    out_M = np.zeros((dim, dim))
    for cnt, i in enumerate(range(dim)):
        out_M[cnt] = M[ind1[cnt]][ind1]
    return(out_M)
#    return(M[ind1][ind2])
#    return(np.transpose(M[ind1])[ind2])


def swap_bas(mol, atom_list, inplace = False):
    """Swap elements in mol._bas that correspond to a pople basis
       to follow the gaussian09(16) indexing.
       If inplace=True, update in situ mol._bas
    """
    swap_list = atom_list
    bas = []
    for cnt, iline in enumerate(mol._atom):
        iat, icrd = iline
        ind = np.where(mol._bas[:,0] == cnt)[0]
        if(iat in swap_list):
            perm_id = perm(ind, ATOMNUM[iat])
            for i in perm_id:
                bas.append(mol._bas[i])
        else:
            for i in ind:
                bas.append(mol._bas[i])
    if(inplace):
        mol._bas = bas
    else:
        pass
    return(bas)

