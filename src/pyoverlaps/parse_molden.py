#!/usr/bin/env python3
#!/usr/bin/env python3
import numpy as np
from math import pi, factorial
from cart2sph import trans_cart2sph as c2s
from cart2sph import offcrt

# Molden parser.
ag2bohr = 1.8897259886 

# Angular momentum numbers
ang_mom = {"s": 0,
           "p": 1,
           "d": 2,
           "f": 3 }

# Number of contractions per shell
ang_mom_ct = [1, 3, 6, 10]
ang_mom_sp = [1, 3, 5, 7]

atom_nums  = {"H": 1, "He": 2, "Li": 3, "Be": 4, "Al": 5, "C": 6,\
        "N": 7, "O": 8, "F": 9, "Ne": 10, "P": 15}

def get_num_ct(bas):
    """Get the number of cartesian functions from a bas array"""
    ncrt = 0
    for cnt, ishl in enumerate(bas):
        ncrt += ang_mom_ct[ishl[1]]
    return(ncrt)

def get_num_sp(bas):
    """Get the number of spherical functions from a bas array"""
    nsp = 0
    for cnt, ishl in enumerate(bas):
        nsp += ang_mom_sp[ishl[1]]
    return(nsp)

def gto_sph_norm(n, a):
    """Prefactor to pseudo-normalize cartesian gtos.
       Missing prefactor (4*pi/(2*n + 1))**0.5.
       The actual multiplication occurs when we transform
       Sxyz into Ssph
    """
    norm = (2**(2*n + 3) * factorial(n+1) * (2*a)**(n + 1.5)\
            /(factorial(2*n + 2) * pi**0.5) )**0.5
    return(norm)

def get_cart2sph(mol):
    """Build cartesian to spherical basis function transformation matrix"""
    bas = mol._bas
    nct = get_num_ct(bas)
    nsp = get_num_sp(bas)
    C   = np.zeros((nsp, nct))

    # row, col counters of C
    i_row = 0
    j_col = 0
    for cn1, ishl in enumerate(bas):
        # Number of spherical/cartesian functions per shell
        isp = ang_mom_sp[ishl[1]]
        ict = ang_mom_ct[ishl[1]]
        # Iterate over rows
        for cn2 in range(isp):
            # Ierate over cols
            for cn3 in range(ict):
                offset = offcrt[ishl[1]]
                c2s_idx = offset + cn3 + cn2*ict
#                print("c2s_idx = ", c2s_idx)
                C[i_row + cn2][j_col + cn3] = c2s[c2s_idx]
#                print("C[" + str(i_row + cn2) + "][" + str(j_col + cn3) + "] = ", c2s[offset + cn3 + cn2*ict] )
        j_col += ict
        i_row += isp
    return(C)

def calc_sph_ovlp_matrix(Sxyz, mol1, mol2 = None):
    """Compute spherical atomic overlap matrix from the (unnormalized)\
    cartesian overlap matrix
    """
    C1 = get_cart2sph(mol1)
    if(mol2 != None):
        C2 = get_cart2sph(mol2)
    else:
        C2 = C1

    Ssph = np.matmul(C1, np.matmul(Sxyz, np.transpose(C2)))
    return(Ssph)

class Mol(object):

    def __init__(self, infile):
        """Mol object, for the time being defined by a molden input"""
        # Lines from infile
        with open(infile, "r") as f:
            self.lines = f.read().splitlines()

        # basic attributes
        self.atomcoords = []
        self.ncrt       = 0 # Number of cartesian basis functions

        # Input arrays for integrals
        self._bas = []
        self._atm = []
        self._env = [0.0 for i in range(0,20)]

    def parse_coords(self):
        for cnt1, iline in enumerate(self.lines):
            if(iline[1:6].lower() == "atoms"):
                if("angs" in iline.lower()):
                    fact = ag2bohr
                elif("bohr" in iline.lower()):
                    fact = 1
                else:
                    fact = 1
                for atid, jline in enumerate(self.lines[cnt1 + 1:]):
                    if(("[" in jline) or (len(jline) == 0)):
                        break
                    ixyz  = jline.split()[-3:] # Last 3 columns
                    ixyz  = [fact * float(i) for i in ixyz]
                    try:
                        atnum = atom_nums[jline.split()[0]] 
                    except:
                        atnum = atom_nums[jline.split()[0][0]] 
                    self.atomcoords.append(ixyz)
                    self._atm.append([atnum, 20 + 4*atid, 1,\
                                      20 + 4*atid + 3, 0, 0]) # add "charge"
                    self._env = self._env + ixyz + [0.0] # add "charge"
                break
#        crds = fact * np.array(crds).astype("float64")
        self.atomcoords = np.array(self.atomcoords).astype("float64")
#        return(crds)

    def parse_shells(self):
        """
        bas: [[atom_idx, amom, n_prm, n_contr, 0, idx_frst_exp,
               idx_first_coeff, 0]..] -> ith shell
        atm: [[nucl_chg, id_x_crd, id+1 z_crd, 0, 0]]
        env: [0, ..., 0, ->0-19 idx
             x_crd_at1, y, z, 0.0, --> x, y, z , buff (chg?) atom 1
             x_crd_at2, y, z, 0.0,
             ...
             1_exp_1_shell, ..., last_exp_1_shell,
             1_cf_1_shell, ..., last_cf_1_shell,
             ...
             1_exp_last_shell, last_exp_last_shell,
             1_cf_last_shell, last_exp_last_shell]
    
        If the same primitives are used for different
        atoms, end saves the primitives only once
        and bas and atm point to these same positions
        """
        prm       = [] # [[exp1, ...,expn, cf1, ..., cfn], ...]
        iatom     = 0 # Atom index
        iprim     = 0 # Primitive index
        skiplines = 0 # Skip primitive lines after reading them 
    #    offset    = 0 # Offset of a given contraction within a shell
        for cnt1, iline in enumerate(self.lines):
            # Pre-processing: read all basis functions
            if(iline[1:4].lower() == "gto"):
                for cnt2, jline in enumerate(self.lines[cnt1 + 1:]):
                    row = jline.split()
                    if(("[" in jline)):
                        break
                    elif(skiplines > 0):
                        skiplines -= 1
                        pass
                    elif(len(row) == 0):
                        pass
                    elif((len(row) >= 1) and (row[0].isdigit())):
                        iatom = int(row[0]) - 1
                    elif(row[0].lower() in ang_mom.keys()):
                        # Construct reference arrays for ith sell
                        Nprim = int(row[1])
                        amom  =  ang_mom[row[0].lower()]
                        self.ncrt  += ang_mom_ct[amom]
                        # Parse all primitives. sp case should also be considered here
                        currid = cnt1 + 1 + cnt2 + 1
                        cf     = []
                        expn   = []
                        for cnt3, kline in enumerate(self.lines[currid: currid + Nprim]):
                            p1, c1 = kline.replace("D","E").split()
                            p1, c1 = float(p1), float(c1)
                            expn.append(p1)
                            # Multiply gto_sph_norm by coeff
                            cf.append(gto_sph_norm(amom,p1)*c1)
#                            cf1.append(c1)
                        if(expn + cf not in prm):
                            prm.append(expn + cf)
                            self._env = self._env + expn + cf
                            idexp     = len(self._env) - 2*Nprim
                            idcf      = len(self._env) - Nprim
                            self._bas.append([iatom, amom, Nprim, 
                                              1, 0, idexp, idcf, 0])
                        else:
                            natom = len(self.atomcoords)
                            id_prm = prm.index(expn + cf)
                            flat_prm = [i for sublist in prm[:id_prm] for i in sublist]
                            idexp = 20 + 4*natom + len(flat_prm)
                            idcf  = idexp + Nprim
#                            idexp = np.where(np.array(self._env) == expn[0])[0][0]
#                            idcf  = idexp + Nprim
                            self._bas.append([iatom, amom, Nprim, 
                                              1, 0, idexp, idcf, 0])
                        skiplines = Nprim
                    
                    elif(row[0].lower() == "sp"): # Consider sp case
                        Nprim = int(row[1])
                        self.ncrt += 4 # + 1s + 3p
                        currid = cnt1 + 1 + cnt2 + 1
                        cf1    = []
                        cf2    = []
                        expn   = []
                        for cnt3, kline in enumerate(self.lines[currid: currid + Nprim]):
                            p1, c1, c2 = kline.replace("D","E").split()
                            p1, c1, c2 = float(p1), float(c1), float(c2)
                            expn.append(p1)
                            # Multiply gto_sph_norm by coeff
                            cf1.append(gto_sph_norm(amom,p1)*c1)
                            cf2.append(gto_sph_norm(amom,p1)*c2)
#                            cf1.append(c1)
#                            cf2.append(c2)
                        iprm = expn + cf1 + expn + cf2 # expn and cf of s and p
                        if(iprm not in prm):
                            prm.append(iprm)
                            self._env = self._env + iprm
                            # Add s shell
                            amom =  ang_mom["s"]
                            idexp     = len(self._env) - 4*Nprim
                            idcf      = len(self._env) - 3*Nprim
                            self._bas.append([iatom, amom, Nprim, 
                                              1, 0, idexp, idcf, 0])
                            # Add p shell
                            amom =  ang_mom["p"]
                            idexp     = len(self._env) - 2*Nprim
                            idcf      = len(self._env) - Nprim
                            self._bas.append([iatom, amom, Nprim, 
                                              1, 0, idexp, idcf, 0])
                        else:
                            # Add s shell
                            amom =  ang_mom["s"]
                            idexp = np.where(np.array(self._env) == iprm[0])[0][0]
                            idcf  = idexp + Nprim
                            self._bas.append([iatom, amom, Nprim, 
                                              1, 0, idexp, idcf, 0])
                            # Add p shell
                            amom =  ang_mom["p"]
                            idexp = idexp + 2*Nprim
                            idcf  = idexp + 3*Nprim
                            self._bas.append([iatom, amom, Nprim, 
                                              1, 0, idexp, idcf, 0])
                        skiplines = Nprim
    
    # MO- referred methods
    def allocate_mo(self, unrestr = False, basis = "spherical"):
        if(basis == "spherical"):
            dim = get_num_sp(self._bas)
        elif(basis == "cartesian"):
            dim = get_num_ct(self._bas)
        else:
            raise ValueError("basis can only be either spherical or cartesian")
        if(unrestr == True):
            self.C_alpha = np.zeros((dim,dim))
            self.C_beta = np.zeros((dim,dim))
            self.C_mo = {"alpha": self.C_alpha,
                      "beta": self.C_beta}
        else:
            self.C_mo = np.zeros((dim,dim))
    
    def get_mo_prop(self):
        """Get MO properties such as Energies,
           spin and occupation numbers
        """
        self.moenergies = []
        self.mooccnos   = [] # MO occupation numbers
        self.mospin     = []
        for cnt, iline in enumerate(self.lines):
            row = iline.split()
            if("Ene=" in iline):
                self.moenergies.append(float(row[1]))
            elif("Spin=" in iline):
                self.mospin.append(row[1])
            elif("Occup=" in iline):
                self.mooccnos.append(float(row[1]))
            else:
                pass
        self.moenergies = np.array(self.moenergies)
        self.mooccnos   = np.array(self.mooccnos)
        self.mospin     = np.array(self.mospin)
    
    def get_all_mo(self, unrestr = False):
        """For now, assume that N_sph has been parsed. In case of unrestricted solutions,
        this method initializes two matrices: C_aplha and C_beta."""
        if(unrestr == True):
            self._get_mo_unr(parser = "Spin= Alpha")
            self._get_mo_unr(parser = "Spin= Beta")
        else:
            self._get_mo()

    def _get_mo(self, parser = "Spin= Alpha"):
        # Parse
        row_cnt = 0
        for cnt1, lin1 in enumerate(self.lines):
#            if(("MO" in lin1)): # or (lin1 == "[MO]")):
            if(("MO" in lin1) and ("MOLDEN" not in lin1)): # or (lin1 == "[MO]")):
                for cnt2, lin2 in enumerate(self.lines[cnt1:]):
                    if("SCF" in lin2): break
                    if(parser in lin2):
                        for cnt3, lin3 in enumerate(self.lines[cnt1 + cnt2 +2:]):
                            if(("Sym=" in lin3) or ( "Ene=" in lin3) or ( "SCF" in lin3) or ("Occup=" in lin3)): break
                            self.C_mo[row_cnt][cnt3] = float(lin3.split()[1])
                        row_cnt += 1

    def _get_mo_unr(self, parser = "Spin= Alpha"):
        spin = parser.replace("Spin= ", "").lower()
        row_cnt = 0
        for cnt1, lin1 in enumerate(self.lines):
#            if(("MO" in lin1)): # or (lin1 == "[MO]")):
            if(("MO" in lin1) and ("MOLDEN" not in lin1)): # or (lin1 == "[MO]")):
                for cnt2, lin2 in enumerate(self.lines[cnt1:]):
                    if("SCF" in lin2): break
                    if(parser in lin2):
                        for cnt3, lin3 in enumerate(self.lines[cnt1 + cnt2 +2:]):
                            if(("Sym=" in lin3) or ( "Ene=" in lin3) or ( "SCF" in lin3)): break
                            self.C_mo[spin][row_cnt][cnt3] = float(lin3.split()[1])
                        row_cnt += 1
    
    def get_mos(self, unrestr = False):
        """Allocate MO arrays and retrieve corresponding values"""
        # Get MOs
        self.allocate_mo(unrestr)
        self.get_all_mo(unrestr)
        #Get properties
        self.get_mo_prop()

    def parse(self, unrestr = False):
        """General parser version"""
        self.parse_coords()
        self.parse_shells()
        self._bas = np.array(self._bas)
        self._env = np.array(self._env)
        self._atm = np.array(self._atm)
#        self.get_mos(unrestr)

if(__name__=="__main__"):
    from argparse import ArgumentParser
    from ovlp_wrapper import calc_cart_ovlp_matrix
    from time import process_time
    from permutations import perm_d_indices, perm_mixed_matrix
    parser = ArgumentParser(description = "Parse molden files. Compute\
            overlap matrix between two sets of molecular orbitals")
    
    parser.add_argument("-f", nargs="+", type = str, \
            help = "Input files")
    opt = parser.parse_args()

    infiles = opt.f

    # Parse molden input file
    mol1 = Mol(infiles[0])
    mol1.parse()
    if(len(infiles)>1):
        mol2 = Mol(infiles[1])
    else:
        mol2 = Mol(infiles[0])
    mol2.parse()
    time1 = process_time()
    #Sxyz
    S    = calc_cart_ovlp_matrix(mol1._bas, mol1._atm, mol1._env, mol1.ncrt,\
                                mol2._bas, mol2._atm, mol2._env, mol2.ncrt) 
    #Ssph
    Ssph = calc_sph_ovlp_matrix(S, mol1, mol2)
    time2 = process_time()
    print("Execution time for cpp_ovlp: ", str(time2 - time1) + " s")

    # Compute MO overlap matrix, if MOs are present
    try:
        mol1.get_mos()
        mol2.get_mos()
        has_mo = True
    except:
        has_mo = False
        pass
    if(has_mo):
        # Permute Ssph so that the ordering of the MOs coincides
        # with that of the molden format
        ind1 = perm_d_indices(mol1)
        ind2 = perm_d_indices(mol2)
        Smat = perm_mixed_matrix(Ssph, ind1, ind2)
        Smo = np.matmul(mol1.C_mo, np.matmul(Smat, np.transpose(mol2.C_mo)))
