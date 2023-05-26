#!/usr/bin/env python3

# Script for aligning two molecular structures
import numpy as np

weights    = {1: 1.0008, 6: 12.011, 7: 14.007, 8: 15.999}
atom_names = {1: "H", 6: "C", 7: "N", 8: "O"}

def quaternion_to_rot(q):
    """Get the 3x3 rotation matrix associated with the 
    input quaternion. q = (qr, qi, qj, qk)"""
    fac = 2./(np.linalg.norm(q)**2)
    print("fac = ", fac)
    d1 = 1 - fac*(q[2]**2 + q[3]**2)
    d2 = 1 - fac*(q[1]**2 + q[3]**2)
    d3 = 1 - fac*(q[1]**2 + q[2]**2)
    # Matrix elements multiplied by the pre-factor
    ri = fac*q[0]*q[1]
    rj = fac*q[0]*q[2]
    rk = fac*q[0]*q[3]
    ij = fac*q[1]*q[2]
    ik = fac*q[1]*q[3]
    jk = fac*q[2]*q[3]

    R = np.array([[d1,    ij - rk, ik + rj],
                  [ij + rk, d2,    jk - ri],
                  [ik - rj, jk + ri, d3   ]])
    return(R)

def translate_matrix(T, Mt):
    """Returns the target geometry (Mt) translated according to to the translation matrix T"""
    return(Mt + T)

def rotate_matrix(Rot, Mt):
    """Returns the target geometry (Mt) rotated according to to the rotation matrix Rot"""
    return(np.matmul(Mt, Rot))

class Align(object):
    """Align the Mtg target matrix to the Mref reference matrix"""
    def __init__(self, Mref, Mtg, atom_list):
        self.A = Mref
        self.B = Mtg
        
        self.dim = len(atom_list)
        self.atom_list = atom_list
        self.atom_weights = [weights[i] for i in self.atom_list]
        self.W = np.zeros((self.dim, self.dim))
        self.K = np.zeros((4,4))

        self.init_matrices()

    def init_matrices(self):
        self.translate()
        self.get_W_matrix()
        self.get_M_matrix()
        self.get_K_matrix()
        self.get_optimal_rotation()

    def align(self):
        """It returns the target coordinates aligned to the reference coordinates"""
#        self.translate()
#        self.get_W_matrix()
#        self.get_M_matrix()
#        self.get_K_matrix()
#        self.get_optimal_rotation()
        R = self.rotate_target()
        return(R)
    
    # Displace COMs
    def translate(self):
        print("Translate called") 
#        CD_A = np.average(self.A, axis = 0, weights = self.atom_weights)
#        print(CD_A)
#        CD_B = np.average(self.B, axis = 0, weights = self.atom_weights)
        CD_A = np.average(self.A, axis = 0)
        print(CD_A)
        CD_B = np.average(self.B, axis = 0)
        self.displ_vec = CD_A - CD_B
        print("displ vec")
        print(self.displ_vec)
        displ = []
        for cnt, i in enumerate(self.B):
            displ.append(self.displ_vec) 
        self.displ = np.array(displ)
        self.B = self.B + self.displ 

#    def translate(self):
#        #print("Translate called")
#   
#        CD_A = (1./len(self.A))*(np.sum(self.A, axis = 0))
#        CD_B = (1./len(self.B))*(np.sum(self.B, axis = 0))
#        self.displ_vec = CD_A - CD_B
#        displ = []
#        for cnt, i in enumerate(self.B):
#            displ.append(self.displ_vec) 
#        self.displ = np.array(displ)
#        self.B = self.B + self.displ

    def get_translation(self):
        """Return displacement vector for target matrix"""
        return(self.displ_vec)

    def get_W_matrix(self):
        """Construct weights matrix. Arg contains a list with the atomic numbers
        of each atom"""
        self.W = np.identity(self.dim)
        for cnt, iat in enumerate(self.atom_list):
            self.W[cnt][cnt] = weights[iat]
    
    def get_M_matrix(self):
        """Construct scalar product matrix:
            M = (A+ * W * B))"""
        self.M = np.matmul(np.transpose(self.A), np.matmul(self.W, self.B))

    def get_K_matrix(self):
        """Construct the matrix K, whose highest eigenvalue determines the 
        optimal rotation matrix"""
        print("getting K")

        # Define explicitly
        
        diag = self.M.diagonal()
        d1 = diag[0] + diag[1] + diag[2]
        d2 = diag[0] - diag[1] - diag[2]
        d3 = -diag[0] + diag[1] - diag[2]
        d4 = -diag[0] - diag[1] + diag[2]
        M = self.M
        # DOUBTS ABOUT THE SIGN OF K[2][1]
        self.K = np.array([[d1, M[1][2] - M[2][1], M[2][0] - M[0][2], M[0][1] - M[1][0]],
                      [M[1][2] - M[2][1], d2, M[0][1] + M[1][0], M[2][0] + M[0][2]],
                      [M[2][0] - M[0][2], M[0][1] + M[1][0], d3, M[1][2] + M[2][1]],
                      [M[0][1] - M[1][0], M[2][0] + M[0][2], M[1][2] + M[2][1], d4]])
        
#        # The following won't work since the matrix K is NOT symmetric!!
##        K = np.zeros((4,4))
##        # Define first row
#        self.K[0][0] = np.trace(self.M)
#        self.K[0][1] = self.M[1][2] - self.M[2][1]
#        self.K[0][2] = self.M[2][0] - self.M[0][2]
#        self.K[0][3] = self.M[0][1] - self.M[1][0]
#
##        # Define its symmetric counterparts
##        K[1][0] = K[0][1]
##        K[2][0] = K[0][2]
##        K[3][0] = K[0][3]
##
##
##        # Define sign permutation matrices
#        P1 = np.array([[0., 0., 1.],
#                       [1., 0., 0.],
#                       [0., 1., 0.]])
#        signs = np.array([-1., -1., 1. ])
#
#        for i in range(1,4):
#            signs = np.matmul(P1,np.transpose(signs))
#            self.K[i][i] = np.matmul(signs, self.M.diagonal())
#            for j in range(i + 1,4):
#                self.K[i][j] = self.M[i-1][j-1] + self.M[j-1][i-1]
##                K[j][i] = K[i][j]
#        self.K = self.K + np.transpose(self.K) - np.matmul(self.K.diagonal(),np.identity(4))
##        
##        
##        print(K)
##        print(" ")
##        print(self.K)

    def get_optimal_rotation(self):
        """Initializes the optimal rotation matrix R"""
        # w = eigenvalues, v = eigenvectors
        w, v = np.linalg.eig(self.K) 
        ind_max = int(np.where(w == w.max())[0])
        # The column selected corresponds to the eigenvector
        # of v.max()
        self.R_arr = v[:,ind_max].astype("float64") 
        self.R = quaternion_to_rot(self.R_arr)
        # In case one uses it externally
#        return(self.R)
    
    def get_rotation_matrix(self):
        """Return optimal rotation matrix"""
        return(self.R)

    def rotate_target(self):
        """Returns the target geometry aligned to the reference geometry"""
#        return(np.transpose(np.matmul(self.R, np.transpose(self.B))))
        return((np.matmul(self.B, self.R)))


#if __name__=="__main__":
#    from parser8 import *
#    from argparse import ArgumentParser
#    parser = ArgumentParser("Read two molden files: a reference and a target,\
#                             then align the target to the reference and print the \
#                             rotated coordinates\
#                            ")
#    parser.add_argument("-f1", dest = "ref", type = str, 
#            help = "Reference geometry file (.molden)")
#    parser.add_argument("-f2", dest = "target", type = str, 
#            help = "Target geometry file (.molden)")
#    options = parser.parse_arguments()
#    infile1 = parser.ref
#    infile2 = parser.target
#
#    # Retrieve geometries
#    mol1 = Parser(infile1)
#    mol1.get_lines()
#    mol1.get_coords_molden()
#    coord1 = np.array(mol1.coords)
#    
#    mol2 = Parser(infile2)
#    mol2.get_lines()
#    mol2.get_coords_molden()
#    coord2 = np.array(mol2.coords)
#
#    atom_list = mol2.atom_list
#
#    # Get rotated geometry
#    rotor = Align(coord1, coord2, atom_list)
#    coord2_new = rotor.align()
#
#    # Print new geometry
##    with open("rotated_geom.molden", "w") as f:
#    with open("rotated_geom.xyz", "w") as f:
#        if(mol2.units == "bohr"):
#            fac = (1./ag2bohr)
#        else:
#            fac = 1.0
##        fac = 1.0
#        f.write(str(mol2.N_atoms) + "\n")
#        f.write("\n")
#        fmt = "{:6s}{:12.6f}{:12.6f}{:12.6f}"
#        for cnt, iat in enumerate(coord2_new):
#            f.write(fmt.format(atom_names[atom_list[cnt]], *(fac*iat)) + "\n")




