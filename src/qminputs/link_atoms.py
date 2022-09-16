# coding: utf-8
from constants import ATOMNUM
import numpy as np
import pytraj as pt
from threading import Thread

def get_unit_vector(ini, fin):
    """Get unit vector along the fin - ini segment, where
    fin = final point (arrow head)
    ini = initial point
    """
    dist = fin - ini
    v =  dist/np.linalg.norm(dist)
    return(v)

class Link_atoms(object):
    def __init__(self, traj, qmmask, igeom = 0, link_atom = "H", link_dist = 1.09):
        """Class to handle eventual link atoms"""

        # GENERAL VARIABLES
#        self.traj       = pt.load(crd, top, frame_indices = [igeom])
        self.traj           = traj
        self.igeom          = igeom
        self.qmmask         = qmmask
        self.qmcrd          = self.traj[qmmask].xyz[igeom]
        self.N_qm           = len(self.qmcrd)
        self.N_mm           = len(self.traj.xyz[igeom]) - self.N_qm
        self.N_link         = 0 # Number of link atoms
        self.link_atom      = link_atom
        self.link_dist      = link_dist # Angstroms
        self.link_atom_def  = "H"  # Default
        self.link_dist_def  = 1.09 # Default

        # INDEXING-RELATED VARIABLES

        self.labels     = {} # {atom_idx: 0 or 1}, 0 = QM atom
                             #                     1 = MM atom
        
        self.arr_bonds  = [] # Nx4 matrix;   [[at1, at2, lab1, lab2],..]
                             # Labels as one-hot encoding defined by the 
                             # labels dictionary
        
        self.qmmm_bonds = [] # Bonds at the QM/MM boundary, with 2 cols:
                             # [[qm_at1, mm_at1],
                             #  [qm_at2, mm_at2],...]
        
        self.nn_atoms   = {} # Dictionary containing all MM nearest 
                             # neighbors (NNs) of each MM atom in a 
                             # QM/MM bond
        self.all_nn     = [] # A list containing all nearest neighbors of all
                             # MM atoms in a QM/MM bond
        
        self.link_pos   = [] # Array containing the positions of all
                             # link atoms
        self.all_link   = [] # All MM atoms in QM/MM boundary + NNs

    def get_link_info(self):
        """Link atom handler"""
        self.get_qm_mm_indices()
        self.get_link_atoms()

    def get_qm_mm_indices(self):
        """Get indices of QM and of MM atoms"""
        
        # Initialize self.labels
        for cnt, iat in enumerate(self.traj.top.atoms):
            # Not so pythonic array comparison
            comp_arr = [(self.traj.xyz[self.igeom][cnt] == iarr).all()\
                        for iarr in list(self.qmcrd)]
            comp     = True in comp_arr
            if(comp):
                self.labels[iat.index] = 0
            else:
                self.labels[iat.index] = 1
#            # HERE'S A BUG!
#            if(self.traj.xyz[self.igeom][cnt] in self.qmcrd):
#                self.labels[iat.index] = 0
#            else:
#                self.labels[iat.index] = 1

        # Initialize self.arr_bonds
        arr_bonds = []
        for i, ib in enumerate(self.traj.top.bonds):
            at0, at1 = ib.indices
            arr_bonds.append([at0, at1,\
                              self.labels[at0], self.labels[at1]])
            
        self.arr_bonds  = np.array(arr_bonds)

    def get_link_atoms(self):
        """Initialize:
            qmmm_bonds: Array with qm/mm bonds
            nn_atoms  : Dictionary with all MM nearest neighbors of MM
                        atoms at the boundary
        """

        # Get atoms involved in QM/MM bonds:
        # Condition for link atoms and indices of QM/MM bonds
        # in arr_bonds
        link_cond  = (((self.arr_bonds[:,2] == 0) &\
                       (self.arr_bonds[:,3] == 1)) |\
                      ((self.arr_bonds[:,2] == 1) &\
                       (self.arr_bonds[:,3] == 0)))
        link_idx   = np.array(np.where(link_cond)).flatten()
        print("link_idx = ", link_idx)
        
        # Iterate over QM/MM bonds only
        for i, ib in enumerate(self.arr_bonds[link_idx]):
            qm_id = int(np.where(ib[2:] == 0)[0])
            mm_id = int(np.where(ib[2:] == 1)[0])
            qm_at = ib[qm_id]
            mm_at = ib[mm_id]
            self.qmmm_bonds.append([qm_at, mm_at])

            # Set link atom positions:
            # Get unit vector along each QM/MM bond
            ini      = self.traj.xyz[self.igeom][qm_at]
            fin      = self.traj.xyz[self.igeom][mm_at]
            v        = get_unit_vector(ini, fin)
            link_crd = ini + self.link_dist * v
            self.link_pos.append(link_crd) 
            
            # Retrieve nearest neighbors of mm_at
            cache_nn = self._get_nn_per_atom(qm_at, mm_at, qm_id)
            self.nn_atoms[mm_at] = cache_nn
            self.all_nn.append(cache_nn)

        self.all_nn = np.array(self.all_nn).flatten()
        self.qmmm_bonds = np.array(self.qmmm_bonds) 
        
        if(len(self.qmmm_bonds)>0):
            # Get number of link atoms
            self.N_link   = len(self.qmmm_bonds)

            # Get all MM atoms involved in QM/MM boundary + NNs
            self.all_link = np.concatenate((self.qmmm_bonds[:,1],\
                                            self.all_nn))
            # Print info regarding the nature of the boundary
            self.print_warning_boundary()
            self.print_warning_charges()
        else:
            print("No link atoms found")

    def print_warning_charges(self):
        """Print a warning message in case the total MM 
           charge of the MM (and QM) region is not an 
           integer
        """
        total_qm = np.sum(self.traj[self.qmmask].top.charge)
        total_mm = np.sum(self.traj["!(" + self.qmmask + ")"].top.charge)

        print("\nWARNING: The total MM charges of the QM region ({:6.4f}) and of the MM region ({:6.4f}) are not integers. Please assess whether to (manually) modify the MM charges so as to obtain an integer value\n".format(total_qm, total_mm))


    def print_warning_boundary(self):
        """Print a warning regarding the link atom and
           the link atom distance in case a non C-C 
           bond is being cut
        """

        # Iterate over qmmm bonds
        for ibond in self.qmmm_bonds:
            at1 = self.traj.top.atom(ibond[0]).atomic_number
            at2 = self.traj.top.atom(ibond[1]).atomic_number
            atname1 = ATOMNUM[at1]
            atname2 = ATOMNUM[at2]
            if((at1 != 6) or (at2 != 6)):
                print("\nWARNING: Cutting a 'non-standard' [C(sp3)-C(sp3)] boundary: The atoms involved are " + atname1 + "(QM) - " + atname2 + "(MM)" )
                if((self.link_atom == self.link_atom_def) and (self.link_dist == self.link_dist_def)):
                    print("Using the DEFAULT parameters:")
                else:
                    print("Using the NON-DEFAULT parameters:")
                print("Link atom:     " + self.link_atom)
                print("Link distance: " + str(self.link_dist))
            else:
                print("\nCutting a 'standard' [C(sp3)-C(sp3)] boundary. The parameters used are:" )
                print("Link atom:     " + self.link_atom)
                print("Link distance: " + str(self.link_dist))

    def _get_nn_per_atom(self, qm_at, mm_at, qm_id):
        """Get nearest neighbors (NN) of mmatom in the MM region
           (that is, all NN of mmatom but qmatom).
           qm_at = QM atom index
           mm_at = MM atom index
        """

        # Find nearest neighbors of MM, excluding the QM atom involved
        # in the QM/MM bond
        nn_cond  = (((self.arr_bonds[:,0] == mm_at) |\
                   (self.arr_bonds[:,1] == mm_at))  &\
                   (self.arr_bonds[:,qm_id] != qm_at))
        
        # Idx of bonds involving mm_at
        nn_id   = np.where(nn_cond)[0] 
        
        cache_nn = []
        for ibond in nn_id:
            i_nn = np.where(self.arr_bonds[ibond][:2] != mm_at)[0]
            cache_nn.append(int(self.arr_bonds[ibond][i_nn]))
        
        return(cache_nn)

    def write_geometry(self, coords_only = True):
        """Write QM and MM coordinates accounting for link atoms. coords_only: whether to print .xyz file for MM positions and charges or just the MM positions"""
        self.write_qm_coordinates()
        self.write_mm_coordinates(coords_only)

    def write_qm_coordinates(self):
        """Write coordinates of atoms in the QM region.\n
        """
        
        # Write QM Coordinates:
        with open("QM_geom{}.xyz".format(self.igeom), "w") as f:
            fmt = "{:<10s}{:>15.7f}{:>15.7f}{:>15.7f}\n"
            qm_atoms = self.N_qm + self.N_link
            f.write("{:<10d}\n\n".format(qm_atoms))
            for nat, iat in enumerate(self.traj.top.atoms):
                if(self.labels[nat] == 0):
                    crd  = self.traj.xyz[self.igeom][nat]
                    name = ATOMNUM[iat.atomic_number]
                    f.write(fmt.format(name, *crd))
                else:
                    pass
            # Print link atoms
            for cnt, icrd in enumerate(self.link_pos):
                f.write(fmt.format(self.link_atom, *icrd))


    def write_mm_coordinates(self, coords_only = True):
        """Write MM point charges and coordinates. Output:\n\
           charges_geomN.xyz = file bearing point charge positions\n\
                               and charges\n\
           MM_geomN.xyz      = file bearing MM positions with atom\n 
                               names (for visualization purposes)\n
        """
        self._write_mm_coordinates()
        if(coords_only):
            pass
        else:
            self._write_mm_charges()

#        try:
#            print("Enabling multithreading to print the MM charges")
#            thread1 = Thread(target = self._write_mm_charges)
#            thread2 = Thread(target = self._write_mm_coordinates)
#            thread1.start()
#            thread2.start()
#            thread1.join()
#            thread2.join()
#        except:
#            print("Multithreading could not be enabled. Printing serially")
#            self._write_mm_charges()
#            self._write_mm_coordinates()


    def _write_mm_charges(self):
        """Write MM point charges"""
        # Write MM point charges
        with open("charges_geom{}.xyz".format(self.igeom), "w") as f:
            fmt = "{:>15.7f}{:>15.7f}{:>15.7f}{:>15.7f}\n"
            MM_atoms = self.N_mm - self.N_link
            f.write("{:<10d}\n\n".format(MM_atoms))
            for nat, iat in enumerate(self.traj.top.atoms):
                if(self.labels[nat] == 1):
                    # Consider only MM atoms
                    if(nat not in self.all_link):
                        # MM atoms not involved in QM/MM boundary
                        crd = self.traj.xyz[self.igeom][nat]
                        chg = self.traj.top.charge[nat]
                        f.write(fmt.format(*crd, chg))
                    elif(nat in self.qmmm_bonds[:,1]):
                        # Print NNs for each QM/MM atom with smeared charges
                        chg_parz = self.traj.top.charge[nat]/\
                                   len(self.nn_atoms[nat])
                        
                        for inn in self.nn_atoms[nat]:
                            # Ith nearest neighbor
                            crd = self.traj.xyz[igeom][inn]
                            chg = self.traj.top.charge[inn] + chg_parz
                            f.write(fmt.format(*crd, chg))
                    else:
                        pass
    
    def _write_mm_coordinates(self):
        """Write M coordinates on an xyz file for visualization
           purposes
        """
        with open("MM_geom{}.xyz".format(self.igeom), "w") as f:
            fmt = "{:<10s}{:>15.7f}{:>15.7f}{:>15.7f}\n"
            MM_atoms = self.N_mm - self.N_link
            f.write("{:<10d}\n\n".format(MM_atoms))
            for nat, iat in enumerate(self.traj.top.atoms):
                if(self.labels[nat] == 1):
                    # Consider only MM atoms
                    if(nat not in self.all_link):
                        # MM atoms not involved in QM/MM boundary
                        crd = self.traj.xyz[self.igeom][nat]
                        name = ATOMNUM[iat.atomic_number]
                        f.write(fmt.format(name, *crd))
                    elif(nat in self.qmmm_bonds[:,1]):
                        # Print NNs for each QM/MM atom
                        
                        for inn in self.nn_atoms[nat]:
                            # Ith nearest neighbor
                            crd    = self.traj.xyz[self.igeom][inn]
                            at_num = self.traj.top.atom(inn).atomic_number
                            name   = ATOMNUM[at_num]
                            f.write(fmt.format(name, *crd))
                    else:
                        pass
    

if(__name__ == "__main__"):
    from argparse import ArgumentParser
    
    # GENERAL ARGUMENTS
    parser = ArgumentParser("Generate xyz coordinates \
            accounting for link atoms.")
    parser.set_defaults(igeom = 0)
    parser.add_argument("-c", dest = "trajectory", type = str,\
                        help = "Trajectory file")
    parser.add_argument("-p", dest = "topology", type = str,\
                        help = "Topology file")
    parser.add_argument("-qm", dest = "qmmask", type = str,\
                        help = "QM Mask")
    parser.add_argument("-ig", dest = "igeom", type = int,\
                        help = "Geometry index")

    # Define variables
    options    = parser.parse_args()
    crd        = options.trajectory
    top        = options.topology
    qmmask     = options.qmmask
    igeom      = int(options.igeom)
    
    traj       = pt.load(crd, top)
#    link_at = Link_atoms(crd, top, qmmask, igeom)
    link_at = Link_atoms(traj, qmmask, igeom)
    link_at.get_link_info()
    link_at.write_geometry(True)
