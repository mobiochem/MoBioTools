#!/usr/bin/env python3 
# Filename: pyoverlaps.py
from align import quaternion_to_rot, translate_matrix, rotate_matrix,\
        Align
import pytraj as pt
from copy import deepcopy
from collections import OrderedDict
from permutations import readcrd, perm, get_d_indices,\
        perm_d_indices, perm_matrix, perm_mixed_matrix, swap_bas
import numpy as np
#import sys

# Global Functions
def get_reference_geometry(infile):
    """Get xyz coordinates of the reference geometry 
    from a molden file. Returns a list of the
    n_qm atomic numbers and a Nx3 array"""
    mol = Parser(infile)
    mol.get_lines()
    mol.get_coords_molden()
    if(mol.units == "bohr"):
        coord = (1./ag2bohr) * np.array(mol.coords) # Transform into Angstrom
    else:
        coord = np.array(mol.coords) # Transform into Angstrom
    return(mol.atom_list, coord)


class Parse_main(object):

    def __init__(self, infile):
        """Class for retrieveing input info."""
        self.main_nsp = {"inp": "", 
                         "ref": "",
                         "crd": "",
                         "top": "",
                         "igeom": 1}
        self.parse(infile)
        

    def parse(self, infile):
        """Parse main input file; update main namespace dictionary"""

        with open(infile, "r") as f:
            lines = f.read().splitlines()
        for line in lines[lines.index("&main") + 1 :]:
            if("&end" in line):
                break;
            row = line.split("=")
            key = row[0].strip()
            val = row[1].strip()
            self.main_nsp[key] = val

class Parse_cas(object):

    def __init__(self, infile):
        """Parse CASSCF input file and create differet 
        namespace disctionaries (only MOLCAS compatible)"""
        # GATEWAY Namespace
        keys_gw = ["title", "coord", "group", "basis set", "xfield"]
        vals_gw = ["default title", "", "C1", "", ""]
        self.gateway_nsp = OrderedDict({i:j for i, j in zip(keys_gw, vals_gw)})
        self.seward_nsp = OrderedDict({})
        self.scf_nsp = OrderedDict({})
        # CASSCF Namespace 
        keys_rs = ["title", "symmetry", "spin", "nactel", "frozen", "inactive", "ras2", "ciroot"]
        vals_rs = ["default title","1", "1", "", "0", "", "", ""]
        self.rasscf_nsp = OrderedDict({i:j for i, j in zip(keys_rs, vals_rs)}) 
        # CASPT2 Namespace
        keys_cp = ["ipea", "shift", "maxiter", "multistate"]
        vals_cp = ["", "", "", ""]
        self.caspt2_nsp = OrderedDict({i:j for i, j in zip(keys_cp, vals_cp)})
        # RASSI Namespace
        self.rassi_nsp = OrderedDict({"nr of jobiph": ""})
        self.parse(infile)

    def parse(self, infile):
        """Parse MOLCAS input file"""
        with open(infile, "r") as f:
            self.lines = f.read().splitlines()

        nspace = [self.lines[cnt].replace("&", "").strip().lower()
                for cnt, line in enumerate(self.lines) if "&" in line]
        print(nspace)

        for name in nspace:
            self.parse_section(name)
    
    def parse_section(self, name):
        """Parse specific section in the molcas namespace;
        e.g: name = gateway"""
        
        uname = "&" + name.upper() # e.g.: &GATEWAY
        print("Parsing " + uname)
        for line in self.lines[(self.lines.index(uname) + 1):]:
            if(("&" in line) or (len(line.split()) == 0)):
                break;
            row = line.split("=")
            self.__getattribute__(name + "_nsp")[row[0].lower().strip()] = row[1].strip()

    def __getattribute__(self, name):
        return(object.__getattribute__(self, name))


class Molcas(object):
    def __init__(self, cas_obj, igeom):
        """Initialize a Molcas object for each geometry to be analyzed.
        Niter = Nth CASSCF iteration after the overlap assessment."""
        self.cas_nsp = deepcopy(cas_obj)
#        self.Niter = Niter
        self.Niter = 1
        self.igeom = igeom
        self.geom = "geom" + str(self.igeom) + ".xyz"
        self.update_cas_nsp()

    def update_cas_nsp(self):
        try:
            self.cas_nsp.gateway_nsp["coord"] = self.geom 
        except:
            raise ValueError("coord not defined in gateway namespace")
        # 14/09/2020: Assessing whether there are point charges
        if(self.cas_nsp.gateway_nsp["xfield"] != ""):
            self.charges = "charges" + str(self.igeom) + ".dat"
            self.cas_nsp.gateway_nsp["xfield"] = self.charges 
        else:
            # 14/09/2020 No point charges considered
            print("No point charges considered")
            del self.cas_nsp.gateway_nsp["xfield"]
#            raise ValueError("xfield not defined in gateway namespace")

    def write_molcas_input(self, alter = "", guessorb = None):
        """Alter is a list of pairs containing the MOs to be swapped"""
        molfile = self.geom.replace(".xyz", "") + "-it" +  str(self.Niter) + ".inp"
        with open(molfile, "w") as f:
            try:
                f.write(">>COPY $CurrDir/" + self.charges + " $WorkDir/." + "\n")
            except:
                pass
            if(self.Niter == 1):
                f.write("&GATEWAY" + "\n")
                for key in self.cas_nsp.gateway_nsp:
                    val = self.cas_nsp.gateway_nsp[key]
                    f.write(key + " = " + val + "\n")
                f.write("RICD" + "\n")
                f.write("&SEWARD &END" + "\n")
                if(guessorb != None):
                    f.write(">>COPY $CurrDir/" + guessorb + " $WorkDir/INPORB"+ "\n")
                else:
                    f.write("&SCF &END" + "\n")
                    f.write(">>COPY $WorkDir/$Project.scf.molden" + " $CurrDir/."+ "\n")
                f.write("&RASSCF &END" + "\n")
                if(guessorb != None):
                    f.write("LumOrb\n")
                for key in self.cas_nsp.rasscf_nsp:
                    val = self.cas_nsp.rasscf_nsp[key]
                    f.write(key + " = " + val + "\n")
                f.write(">>COPY $WorkDir/$Project.rasscf.molden" + " $CurrDir/."+ "\n")
                f.write(">>COPY $WorkDir/$Project.RasOrb" + " $CurrDir/."+ "\n")
                self.Niter +=1
            else:
                f.write("&GATEWAY" + "\n")
                for key in self.cas_nsp.gateway_nsp:
                    val = self.cas_nsp.gateway_nsp[key]
                    f.write(key + " = " + val + "\n")
                f.write("RICD" + "\n")
                f.write("&SEWARD &END" + "\n")
                guessorb= self.geom.replace(".xyz", "") + "-it" + str(self.Niter - 1) + ".inp.RasOrb"
                f.write(">>COPY $CurrDir/" +guessorb + " $WorkDir/INPORB" + "\n")
                f.write("&RASSCF &END" + "\n")
                for key in self.cas_nsp.rasscf_nsp:
                    val = self.cas_nsp.rasscf_nsp[key]
                    f.write(key + " = " + val + "\n")
                f.write("alter" + "\n")
                f.write(str(len(alter)) + "\n")
                for cnt, ipair in  enumerate(alter):
                    f.write("{:<2d} {:<3d} {:<3d}".format(1, *ipair) + "\n")
                f.write("LumOrb" + "\n")

                f.write(">>COPY $WorkDir/$Project.rasscf.molden" + " $CurrDir/."+ "\n")
                f.write(">>COPY $WorkDir/$Project.RasOrb" + " $CurrDir/."+ "\n")
                self.Niter +=1

class Parse_crd(object):

    def __init__(self, ref_coord, atom_list, traj, qmmask, top = None):
        """Parse MD trajectory. 
        ref_coord is a matrix with the reference coordinates.\
        atom_list is a list with the atomic number of each qm atom"""
        self.ref_coord = ref_coord
        self.n_qm = len(self.ref_coord)
        self.atom_list = atom_list
        if(top != None):
            self.traj = pt.load(traj, top)
        self.charges = self.traj.top.charge
        self.n_atoms = self.traj.top.n_atoms
        self.n_charges = self.n_atoms - self.n_qm
        self.qmmask = qmmask

    def rotate_crd(self, igeom):
        """Rotate ith geometry. It returns the ith 
        geometry rotated (Nx3 matrix)"""
        #qm_coord = self.traj[igeom][0:self.n_qm]
        qm_coord = self.traj[self.qmmask][igeom].xyz
        rotor = Align(self.ref_coord, qm_coord, self.atom_list)
        # Translate matrix
        Tn = translate_matrix(rotor.displ_vec, self.traj[igeom][:])
        # Rotate matrix
        Rot = rotor.get_rotation_matrix()
        M = rotate_matrix(Rot, Tn)
        return(M)

    def write_crd(self, igeom):
        """Write xyz coordinates of the qm atoms of the igeom snapshot
        of self.traj. The coordinates need to be first rotated with respect
        to the ref_coord array."""
        # Get rotated matrix
        M = self.rotate_crd(igeom)

        # Get index of residue in qmmask
        for ires in self.traj.top.residues:
            if(ires.name in self.qmmask):
                self.ires = ires.index
                break
        print("ires = ", self.ires)

        # Retrieve atoms in qmmask
        atoms_qmmask = []
        for iat in self.traj.top.atoms:
            if(iat.resid == self.ires):
                atoms_qmmask.append(iat)

        # Write qm geometry
        with open("geom" + str(igeom) + ".xyz", "w") as f:
            f.write(str(self.n_qm) + "\n")
            f.write("\n")
            fmt1 = "{:<7s}{:<15.6f}{:<15.6f}{:<15.6f}"
            for cnt, iat in enumerate(self.traj.top.atoms):
                if(iat.resid == self.ires):
                    f.write(fmt1.format(iat.name, *M[cnt][0:3]) + "\n")
                    
        if(len(self.charges)>0):
            with open("charges" + str(igeom) + ".dat", "w") as f:
                f.write(str(self.n_charges) + "\n")
                fmt2 = "{:<15.6f}{:<15.6f}{:<15.6f}{:<15.6f}{:<10.1f}{:<10.1f}{:<10.1f}"
                for cnt, iat in enumerate(self.traj.top.atoms):
                    if(iat.resid != self.ires):
                        f.write(fmt2.format(M[cnt][0], M[cnt][1], M[cnt][2],\
                            self.charges[cnt], 0.0, 0.0, 0.0) + "\n")

def gen_submission():
    with open("tmp_env_vars.sh", "w") as f:
        txt = """#!/bin/sh
export JOBID=$$
export Project=$1
if [ -z $WorkDir ]; then
 export WorkDir=$MOLCAS_SCRATCH/$USER/$Project-${JOBID}
 mkdir -p $WorkDir
#  echo created temp directory on $WorkDir
 if [ ! -d $WorkDir ]; then
  export WorkDir=$PWD/TEMP
#  echo created temp directory on $WorkDir
  mkdir -p $WorkDir
 fi
fi


CurrDir=$PWD

if [ -z $MOLCASMEM ]; then
 export MOLCASMEM=6000
fi
 export MOLCAS_MEM=$MOLCASMEM
 export MOLCAS_MAXMEM=$MOLCASMEM
if [ -f $Project ]
then
    cp $CurrDir/$Project $WorkDir/.
fi

# Run openmolcas
$runmolcas $Project > $Project.out
if [ -e $WorkDir/$Project.JobIph ]; then
 bzip2 < $WorkDir/$Project.JobIph > $CurrDir/$Project.JobIph.bz2
fi

# Retrieve output data
cd $CurrDir
rm -rf $WorkDir
"""
        f.write(txt)
        comm = "chmod a+rwx " + os.getcwd() + "/tmp_env_vars.sh"
        subprocess.call(comm.split())

if __name__ == "__main__":
    import pickle
    import os
    from pathlib import Path
    import subprocess
    from ovlp_wrapper import symmetrize, wrap_ovlp_matrix, wrap_ovlp_matrix_mix,\
            calc_cart_ovlp_matrix
    from parse_molden import Mol, calc_sph_ovlp_matrix, \
            get_cart2sph, gto_sph_norm, get_num_sp, ag2bohr
    from cart2sph import trans_cart2sph as c2s
    from cart2sph import offcrt
    from align import Align
    from matrixio import save_matrix
    from time import process_time
    from argparse import ArgumentParser

    # Input arguments
    parser = ArgumentParser("PyOverlaps MO recovery pof the active space")
    parser.set_defaults(units = "bohr", calculate = "yes",\
            restart = 1, tpl = "tpl.inp", ngeom = 1, igeom = 0, align = True, guess=None)
    parser.add_argument("-p ", dest = "top", type = str,\
            help = "Topology file")
    parser.add_argument("-a ", dest = "align", type = bool,\
            help = "Align geometries?")
    parser.add_argument("-c", dest = "crd", type = str,\
            help = "Trajectory file")
    parser.add_argument("-r", dest = "reference", type = str,\
            help = "Reference molden file")
    parser.add_argument("-tpl ", dest = "template", type = str,\
            help = "Molcas template file")
    parser.add_argument("-rng", nargs = "+", type=int,\
            help = "Range of MOs in the active space (e.g: 30 43)")
#    parser.add_argument("-inf", dest = "infimum", type = str,
#            help = "Minimum MO label to consider")
#    parser.add_argument("-sup", dest = "supremum", type = str,
#            help = "Maximum MO label to consider")
    parser.add_argument("-cl", dest = "calc", type = str,\
            help = "Perform CASSCF intermediate calculations? (yes/no)")
    parser.add_argument("-rst", dest = "restart", type = int,\
            help = "Restart from iteration N (N>=1)")
    parser.add_argument("-u", dest = "units", type = str, 
            help = "Units (bohr, angstrom)")
#    parser.add_argument("-b", dest = "basis", type = str,
#            help = "Basis type, either cartesian or spherical")
    parser.add_argument("-qm", dest = "qmmask", type = str,
            help = "QM mask")
    parser.add_argument("-g", dest = "guess", type = str,
            help = "Initial guess of molecular orbitals (usually same as reference)")
    parser.add_argument("--ig", dest = "igeom", type = int, help = "Specific\
                 frame for which to generate an input file. Default = 0")
    parser.add_argument("--n", dest = "ngeom", type = int, help = "Number of\n\
                         geometries to generate. If chosen, do not set --ig. Default = 1")

    options    = parser.parse_args()
    trajfile   = options.crd
    topfile    = options.top
    reference  = options.reference
    tpl        = options.template
    inf, sup   = options.rng[0:2]       
    rst        = options.restart
    calculate  = options.calc
    qmmask     = options.qmmask
    igeom      = options.igeom
    ngeom      = options.ngeom
    align      = options.align
    guess      = options.guess
    
    units      = options.units
    basis      = "spherical"
    links      = [trajfile, topfile, tpl, reference]
    if(guess != None):
        links.append(guess)

    inf       = inf - 1 # Internal indexing
#    sup       = sup - 1
#    active_sp = list(range(inf, sup+1)) # Indexing starts from zero
    active_sp = list(range(inf, sup)) # Indexing starts from zero
    maxit = 5 # Max number of CASSCF iterations
    workdir = os.getcwd() + "/"

    #Create working directory for the given igeom
    idir = workdir + "geom" + str(igeom)
    pth  = Path(idir)
    pth.mkdir(parents = True, exist_ok = True)
    for ilink in links:
        try:
            os.symlink(src = workdir + ilink, dst = idir + "/" + ilink)
        except:
            pass
    os.chdir(idir)
    # Generate folder containing temporary matrices
    comm = "mkdir -p matrices"
    subprocess.call(comm.split())

    # Start parsing
    cas_obj = Parse_cas(tpl)
    molcas_obj = Molcas(cas_obj, igeom)
    molcas_obj.write_molcas_input(guessorb = guess)
    converged = False

    # Parse coordinates:
    # Retrieve reference geometry info
    print(idir + "/" + reference)
    mol1      = Mol(idir + "/" + reference)
    mol1.parse()
    ref_coord =  mol1.atomcoords
    atom_list = mol1._atm[:,0]
    # This performs automatically the alignment
    if(units == "bohr"):
        ref_coord = (1/ag2bohr) * ref_coord
    crd_obj   = Parse_crd(ref_coord,
                          atom_list,
                          trajfile,
                          qmmask,
                          topfile) 
    crd_obj.write_crd(igeom)

    
    # Analyze MO overlaps for a given geometry
    print("Entering MO overlap analysis")
    print("Reference geometry taken from " + reference)
    print("Reference active space: " + str(inf + 1) + "-" + str(sup))
    for it in range(rst, maxit + 1):
        mos_to_add = [] # List of MOs to be put in the AS
        idx_add = [] # List of MOs to be put in the AS
        mos_to_del = [] # List of MOs to be removed from the AS
        idx_del = [] # List of MOs to be removed from the AS
        min_lb = 0
#        max_lb = 2 * sup
        max_lb = int(1.5 * sup)
        print("Iteration " + str(it))
        # 1. Call MOLCAS
        inp_cas = "geom" + str(igeom) + "-it" + str(it) + ".inp"
        time1 = process_time()
        if(calculate == "yes"):
            print("Executing MOLCAS: input file " + inp_cas)
            gen_submission()
            comm = os.getcwd() + "/tmp_env_vars.sh " + inp_cas + " 1"
            try:
                subprocess.call("echo Now executing command "+comm, shell = True)
                subprocess.call(comm.split())
            except:
                raise OSError("molcas exit wit error status 1")
            time2 = process_time()
            print("Execution time for MOLCAS: ", str(time2 - time1) + " s")
        else:
            print("Skipping MOLCAS execution")
       
       
       # 2. Calculate MO overlap matrix
        print("MO analysis for iteration " + str(it))

        mol2 = Mol(inp_cas + ".rasscf.molden")
        mol2.parse()
        # Align snapshot coordinates to those of the reference

        # Calculate atomic overlap matrix
        Sxyz   = calc_cart_ovlp_matrix(mol1._bas, mol1._atm, mol1._env, mol1.ncrt,\
                                       mol2._bas, mol2._atm, mol2._env, mol2.ncrt) 
        Ssph   = calc_sph_ovlp_matrix(Sxyz, mol1, mol2)
        if(basis == "spherical"):
            # Consider permutation of d orbitals (pyscf to gaussian numbering)
            # Here we are assuming that mol1 and mol2 bases have the same numbering
            ind1 = perm_d_indices(mol1)
            ind2 = perm_d_indices(mol2)
#            Smat = perm_matrix(Ssph, ind)
            Smat = perm_mixed_matrix(Ssph, ind1, ind2)

        else:
            Smat = Sxyz # Warning: Unnormalized

        print("Calculation of MO overlap matrix")
        # Calculating molecular orbital overlap matrix
        mol1.get_mos()
        mol2.get_mos()
        OVL = np.matmul(mol1.C_mo, np.matmul(Smat, np.transpose(mol2.C_mo)))
        OVL_tr = OVL[inf:sup,inf:sup]
        
        # Retrieve indices of the maximum value for each column
        # i.e: for each sample MO
        print("analyzing MOs within the range: ", str(min_lb) + "-" + str(max_lb))
        warning = 0
        corr_mo = {} 

        # 14/09/2020: MOs that could be potentially added because their col maximum coinciced
        # with a reference active space MO
        pot_add_col = []
        pot_add_row = []
        # Iteration over columns within range
        # 02/09/2020 Consider full range of columns
        for colid in range(len(OVL)):
            # Row index of the maximum
            maxid = np.abs(OVL[:,colid]).argmax()
            max_val = np.abs(OVL[:,colid]).max()
            max_as = np.abs(OVL[inf:sup,colid]).max() # Maximum within the AS
            threshold = 0.02 # 04/09/2020
            try:
                prev_alter = alter
                prev_equil = equil
            except:
                prev_alter = []
                prev_equil = []
            print("idx: ", maxid + 1, colid + 1, "val: ", OVL[maxid, colid])
            if(maxid in active_sp):
                # See whether this MO should be added
                if (colid in active_sp):
                    pass
                else:
                    # 08/09/2020: Guarantee uniqueness of the ref MOs
                    # 08/09/2020: Variables borrowed from mos_to_add>mos_to_del criterion 
                    row_max_id = np.abs(OVL[maxid,:]).argmax() # Column label of the row maximum
                    row_max = np.abs(OVL[maxid,:]).max()
                    if(maxid + 1 in idx_add):
                        print("maxid = ", maxid + 1, "in idx_add: add to mos_to_add only the one with the highest overlap")
                        # Add only the one with the highest overlap
                        row_max_pos = np.where(idx_add == maxid + 1)[0][0] # Position on idx_add
                        col_max = mos_to_add[row_max_pos] - 1 # column maximum relative to that position
                        max_rel = max(np.abs(OVL[maxid][col_max]), np.abs(OVL[maxid][colid]))
                        max_rel_col = np.where(np.abs(OVL[maxid,:]) == max_rel)[0][0] # Column index of the maximum between both candidates
                        print("col_max = ", col_max + 1) 
                        print("colid = ", colid + 1)
                        print("max_rel_col = ", max_rel_col + 1)
                        mos_to_add = np.delete(mos_to_add, row_max_pos).tolist()
                        idx_add = np.delete(idx_add, row_max_pos).tolist()
                        mos_to_add.append(max_rel_col + 1)
                        idx_add.append(maxid + 1)
                    # 08/09/2020: Guarantee that the max_val is greater than the row maximum,
                    # and that the row maximum is not in the active space (row criterion).
                    # Same as row criterion for mos_to_add>mos_to_del
                    elif(row_max_id in active_sp):
                        # Do not add the element. 
                        # 14/09/2020 MO to be potentially added
                        pot_add_col.append(colid + 1)
                        pot_add_row.append(maxid + 1)
                        pass
                    else:
                        mos_to_add.append(colid + 1)
                        idx_add.append(maxid + 1)
                # Update corr_mo
                corr_mo[maxid + 1] = colid + 1
            else:
                # See whether this MO should be deleted
                if(colid in active_sp):
                    # 15/09/2020: Guarantee uniqueness of the ref MOs
                    if(maxid + 1 in idx_del):
                        print("maxid = ", maxid + 1, "in idx_del: add to mos_to_del only the one with the highest overlap")
                        # Delete only the one with the highest overlap
                        row_max_pos = np.where(idx_del == maxid + 1)[0][0] # Position on idx_del
                        col_max = mos_to_del[row_max_pos] - 1 # column maximum relative to that position
                        max_rel = max(np.abs(OVL[maxid][col_max]), np.abs(OVL[maxid][colid]))
                        max_rel_col = np.where(np.abs(OVL[maxid,:]) == max_rel)[0][0] # Column index of the maximum between both candidates
                        print("col_max = ", col_max + 1) 
                        print("colid = ", colid + 1)
                        print("max_rel_col = ", max_rel_col + 1)
                        mos_to_del = np.delete(mos_to_del, row_max_pos).tolist()
                        idx_del = np.delete(idx_del, row_max_pos).tolist()
                        mos_to_del.append(max_rel_col + 1)
                        idx_del.append(maxid + 1)
                    elif(abs(max_val - max_as) < threshold):
                        print("Not swapping: " + str(colid + 1) + " max_as = " + "{:10.5f}".format(max_as))
                        # Do not swap if there's a high ovl also within the ref.
                        # Active space
                        pass
                    else:
                        mos_to_del.append(colid + 1)
                        idx_del.append(maxid + 1)
                else:
                    pass
        # End of iteration over columns within range
        print("corr_mo after main algorithm: ", sorted(corr_mo.items()))
        # 04/09/2020: Ensure that all corr_mo are assigned before assessing the lengths
        # of mos_to_add and mos_to_del
        # Assess whether both mos_to_add and mos_to_del have the same length
        missing_mos = np.setdiff1d(list(range(inf + 1, sup + 1)), list(corr_mo.keys()))
        if (len(missing_mos) != 0):
            print("corr_mo does not have all ref MOs assigned. \n An update attempt will be done \n" +  
            "avoiding repetitions. Use max row criterion")
            print("missing_mos = ", missing_mos)
            for imiss in missing_mos:
                colmax = np.abs(OVL[imiss - 1,:]).argmax()
                max2 = np.sort(np.abs(OVL[imiss - 1,:]))[-2] # Second highest value in row
                max3 = np.sort(np.abs(OVL[imiss - 1,:]))[-3] # Second highest value in row
                max2_in = max2 in OVL[imiss - 1,:]
                max3_in = max3 in OVL[imiss - 1,:]
                if(max2_in):
                    col2max = np.where(OVL[imiss - 1,:] == max2)[0][0]
                else:
                    col2max = np.where(OVL[imiss - 1,:] == -max2)[0][0]
                if(max3_in):
                    col3max = np.where(OVL[imiss - 1,:] == max3)[0][0]
                else:
                    col3max = np.where(OVL[imiss - 1,:] == -max3)[0][0]

                # Assess whether updating mos_to_add/mos_to_del
                if(colmax + 1 not in list(corr_mo.values())):
                    print("assessing colmax = ", colmax + 1)
                    if((colmax not in active_sp) and (colmax + 1 not in mos_to_add)):
                        mos_to_add.append(colmax + 1)
                        idx_add.append(imiss)
                        corr_mo[imiss] = colmax
                        print("imiss=", imiss, " colmax = ", colmax + 1)
                    elif((colmax in active_sp) and (colmax + 1 not in mos_to_del)):
                        # If it in mos_to_del do nothing, the assessment occurs below (since it is more complicated)
                        corr_mo[imiss] = colmax
                        print("imiss=", imiss, " colmax = ", colmax + 1)

                elif(col2max + 1 not in list(corr_mo.values())):
                    print("assessing col2max = ", col2max + 1)
                    if((col2max not in active_sp) and (col2max + 1 not in mos_to_add)):
                        mos_to_add.append(col2max + 1)
                        idx_add.append(imiss)
                        corr_mo[imiss] = col2max
                        print("imiss=", imiss, " col2max = ", col2max + 1)
                    elif((col2max in active_sp) and (col2max + 1 not in mos_to_del)):
                        # If it in mos_to_del do nothing, the assessment occurs below (since it is more complicated)
                        corr_mo[imiss] = col2max
                        print("imiss=", imiss, " col2max = ", col2max + 1)

                elif(col3max + 1 not in list(corr_mo.values())):
                    print("assessing col3max = ", col3max + 1)
                    if((col3max not in active_sp) and (col3max + 1 not in mos_to_add)):
                        mos_to_add.append(col3max + 1)
                        idx_add.append(imiss)
                        corr_mo[imiss] = col3max
                        print("imiss=", imiss, " col3max = ", col3max + 1)
                    elif((col3max in active_sp) and (col3max + 1 not in mos_to_del)):
                        # If it in mos_to_del do nothing, the assessment occurs below (since it is more complicated)
                        corr_mo[imiss] = col3max
                        print("imiss=", imiss, " col3max = ", col3max + 1)
                else:
                    print("corr_mo not updated\n ")
            print("After corr_mo correction:")
        print("mos to add: ", "row idx: ", idx_add, "col idx: ", str(mos_to_add))
        print("mos to del: ", "row idx: ", idx_del, "col idx: ", str(mos_to_del))
        print()
            
        # 04/09/2020: END Ensure that all corr_mo are assigned before assessing the lengths

        if(len(mos_to_add) != len(mos_to_del)):
#        while(len(mos_to_add) != len(mos_to_del)):
            if(len(mos_to_add) < len(mos_to_del)):
                print("Warning, mos_to_add < mos_to_del; row criterion applied")
#                warning = 1
                # A MO outside of the active space needs to be recovered,
                # but its maximum ovl lies outside the active space of ref.
                # It could potentially have large contributions from MOs within
                # the reference AS.
                # Strategy 1: See whether the corr_mo is lacking a MO from the ref AS,
                # then analyze the maximum overlap along the row and see whether it is
                # outside of the sample active space
                missing_mos = np.setdiff1d(list(range(inf + 1, sup + 1)), list(corr_mo.keys()))
                print("missing_mos = ", missing_mos)
                for imiss in missing_mos:
                    colmax = np.abs(OVL[imiss - 1,:]).argmax()
                    max2 = np.sort(np.abs(OVL[imiss - 1,:]))[-2] # Second highest value in row
                    max3 = np.sort(np.abs(OVL[imiss - 1,:]))[-3] # Second highest value in row
                    max2_in = max2 in OVL[imiss - 1,:]
                    max3_in = max3 in OVL[imiss - 1,:]
                    if(max2_in):
                        col2max = np.where(OVL[imiss - 1,:] == max2)[0][0]
                    else:
                        col2max = np.where(OVL[imiss - 1,:] == -max2)[0][0]
                    if(max3_in):
                        col3max = np.where(OVL[imiss - 1,:] == max3)[0][0]
                    else:
                        col3max = np.where(OVL[imiss - 1,:] == -max3)[0][0]

                    if((colmax not in active_sp) and (colmax + 1 not in mos_to_add)):
                        mos_to_add.append(colmax + 1)
                        idx_add.append(imiss)
                        corr_mo[imiss] = colmax
                        print("imiss=", imiss, " colmax = ", colmax + 1)
                    elif((col2max not in active_sp) and (col2max + 1 not in mos_to_add)):
                        # Attempt to find the column idx using the seconod highest row value
                        mos_to_add.append(col2max + 1)
                        idx_add.append(imiss)
                        corr_mo[imiss] = col2max
                        print("imiss=", imiss, " col2max = ", col2max + 1)
                    elif((col3max not in active_sp) and (col3max + 1 not in mos_to_add) ):
                        # Attempt to find the column idx using the seconod highest row value
                        mos_to_add.append(col3max + 1)
                        idx_add.append(imiss)
                        corr_mo[imiss] = col3max
                        print("imiss=", imiss, " col2max = ", col3max + 1)
                    else:
                        warning = 1

                # 14/09/2020: Consider eventual MOs that were not added at the beginning
                if((len(missing_mos) == 0) or (warning == 1)):
                    for cnt, ipot in enumerate(pot_add_col):
                        if(len(mos_to_add) == len(mos_to_del)):
                            break
                        if(ipot not in mos_to_add):
                            mos_to_add.append(ipot)
                            idx_add.append(pot_add_row[cnt])

            else:
                # A column is equivalent to a reference MO, however it might be that there is another sampled MO
                # within the active space that is also equivalent to that reference MO. This implies that the corr_mo
                # dictionary has all ref MOs assigned. 
                print("Warning, mos_to_add > mos_to_del; Second row criterion applied")
                #warning = 1
                print()
                diff_len = len(mos_to_add) - len(mos_to_del) # >=1
                id_to_del = [] # Eventual indices of elements of mos_to_add that need to be deleted
                for i in range(len(mos_to_del), len(mos_to_add)):
                    row_idx = idx_add[i] # 07/09/2020: Iterate from the first extea element of mos_to_add
                    col_idx = mos_to_add[i]
                    # See whether the row maximum happens to be within the sample active space 
                    col_eq = corr_mo[row_idx] # It need not be within the sample active space!
                    row_max_id = np.abs(OVL[row_idx - 1,:]).argmax()
                    row_max = np.abs(OVL[row_idx - 1,:]).max()
                    if(row_max_id in active_sp):
                        # This implies that the column maximum is >= row_max > ovl. of mo_to_add
                        # and thus, the MO should not be added
                        print("deleting element ", i, " from mos_to_add")
                        id_to_del.append(i) # 07/09/2020: indices of elem. of mos_to_add to be deleted
                        print("id_to_del = ", id_to_del)
                    # The row maximum is outside of the active space. We need to find
                    # a candidate to remove from the active space. First candidate: col_eq
                    elif(col_eq - 1 in active_sp):
                        # This automatically implies that the maximum of that column has ovl< row_max,
                        # and therefore should be deleted
                        mos_to_del.append(col_eq)
                        idx_del.append(row_idx)
                    # We need to find another candidate to remove from the active space:
                    # Find the col within the active space that has the max ovl with row_idx 
                    else:
                        # Its value is clearly < row_max
                        row_max_as = np.abs(OVL[row_idx - 1,inf:sup]).max()
                        row_max_in = row_max_as in OVL[row_idx - 1,:]
                        if(row_max_in):
                            max_as = np.where(OVL[row_idx - 1,:] == row_max_as)[0][0]
                        else:
                            max_as = np.where(OVL[row_idx - 1,:] == -row_max_as)[0][0]
                        mos_to_del.append(max_as + 1)
                        idx_del.append(row_idx)
                        # Eventually see whether it is its column maximum
                # End iteration, now delete eventual mos_to_add elements
                mos_to_add = np.delete(mos_to_add, id_to_del).tolist()
                idx_add = np.delete(idx_add, id_to_del).tolist()




        # END 28/04/2020
        print("After mos_to_add>=<mos_to_del assessment:")
        print("mos to add: ", "row idx: ", idx_add, "col idx: ", str(mos_to_add))
        print("mos to del: ", "row idx: ", idx_del, "col idx: ", str(mos_to_del))
        # 03/09/2020
        alter = list(zip(mos_to_add, mos_to_del)) # Couples of MOs to be swapped
        equil = list(zip(idx_add, idx_del)) # To ensure that we do not count successful false positives
        mo_swap_len = max(len(mos_to_add), len(mos_to_del))
        print("MOs to swap at current iteration:")
        for cnt, ipair in enumerate(alter):
            print("{:3d}{:3d}".format(*ipair))

        # Save AO and MO overlap matrix at every step
        save_matrix(OVL, "matrices/overlaps_mo_full_" + str(it), ".dat", symmetric = False)
        save_matrix(OVL_tr, "matrices/overlaps_mo_tr_" + str(it), ".dat", symmetric = False)
        with open("matrices/OVL_" + str(it) +  ".bin", "wb") as f:
            pickle.dump(OVL, f)
        with open("matrices/Smat_" + str(it) +  ".bin", "wb") as f:
            pickle.dump(Smat, f)

        if((len(alter) == 0) and (mo_swap_len == 0)):
            converged = True
            print("Convergence attained at step " + str(it))
            print("Saving MO overlap matrix in overlaps_mo_full.dat")
            save_matrix(OVL, "overlaps_mo_full", ".dat", symmetric = False)
            save_matrix(OVL_tr, "overlaps_mo_tr", ".dat", symmetric = False)
            with open("OVL.bin", "wb") as f:
                pickle.dump(OVL, f)
            print("Final iteration = " + str(it))
            print("MOs swapped at the previous iteration")
            for cnt, ipair in enumerate(prev_alter):
                print("{:3d}{:3d}".format(*ipair))
            print("Corresponding MOs in the reference set")
            for cnt, ipair in enumerate(prev_equil):
                print("{:3d}{:3d}".format(*ipair))
            # Print summary file, in case N iterations >1
            print("Printing out summary file...")
            with open("summary.dat", "w") as f:
                formtit = "{:<10s}{:<20s}{:<20s}{:<20s}{:<20s}"
                form = "{:<10d}{:<20s}{:<20s}{:<20d}{:<20d}"
                f.write(formtit.format("Geom", "Last MOs swapped", "Corresp. ref. MOs", "N iterations", "Warnings") + "\n")
                print(formtit.format("Geom", "Last MOs swapped", "Corresp. ref. MOs", "N iterations", "Warnings"))
                for cnt in range(len(prev_alter)):
                    alt = str(prev_alter[cnt][0]) + " " + str(prev_alter[cnt][1])
                    ref = str(prev_equil[cnt][0]) + " " + str(prev_equil[cnt][1])
                    f.write(form.format(igeom, alt, ref, it, warning) + "\n")
                    print(form.format(igeom, alt, ref, it, warning))
            break;
        elif((converged == False) and (it == maxit)):
            print("The maximum number of iterations (" + str(maxit) + ") \
                    has been attained \n without achieveing convergence.")
            print("Saving MO overlap matrix in overlaps_mo_full.dat")
            save_matrix(OVL, "overlaps_mo_full", ".dat", symmetric = False)
            save_matrix(OVL_tr, "overlaps_mo_tr", ".dat", symmetric = False)
            with open("OVL.dat", "wb") as f:
                pickle.dump(OVL, f)
            print("Final iteration = " + str(it))
            print("Convergence not attained at step " + str(it))
            print("MOs swapped at the previous iteration")
            for cnt, ipair in enumerate(prev_alter):
                print("{:3d}{:3d}".format(*ipair))
            print("Corresponding MOs in the reference set")
            for cnt, ipair in enumerate(prev_equil):
                print("{:3d}{:3d}".format(*ipair))
            
            # Print summary file
            with open("summary.dat", "w") as f:
                formtit = "{:<10s}{:<20s}{:<20s}{:<20s}{:<20s}"
                form = "{:<10d}{:<20s}{:<20s}{:<20d}{:<20d}"
                f.write(formtit.format("Geom", "Last MOs swapped", "Corresp. ref. MOs", "N iterations", "Warnings") + "\n")
                print(formtit.format("Geom", "Last MOs swapped", "Corresp. ref. MOs", "N iterations", "Warnings"))
                for cnt in range(len(prev_alter)):
                    alt = str(prev_alter[cnt][0]) + " " + str(prev_alter[cnt][1])
                    ref = str(prev_equil[cnt][0]) + " " + str(prev_equil[cnt][1])
                    f.write(form.format(igeom, alt, ref, it, warning) + "\n")
                    print(form.format(igeom, alt, ref, it, warning))
        else:
            molcas_obj.write_molcas_input(alter = alter)
    
    os.chdir(workdir)
