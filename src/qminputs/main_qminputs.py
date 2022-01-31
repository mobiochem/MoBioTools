#!/usr/bin/env python3
import numpy as np
import pytraj as pt
from copy import deepcopy
from collections import OrderedDict
from argparse import ArgumentParser
import os
from common_parser import QMtemplate
from nwchem_parser import NWChemtpl
from gau_parser import Gautpl
from molcas_parser import Molcastpl

class Main(object):

    def __init__(self, infile):
        """Class that parses the main input file"""
        self.namelist = {"tpl": None,
                          "top": None,
                          "traj": None,
                          "geoms": None,
                          "ref": None,
                          "solvmask": None,
                          "closest": None,
                          "bsse": None}
        for key in self.namelist.keys():
            self[key] = self.namelist[key]
        self._infile   = infile
        self.parse()

    def parse(self):
        """Parse main input file"""

        with open(self._infile, "r") as f:
            lines = f.read().splitlines()

        for line in lines[lines.index("&main") + 1 :]:
            if("&end" in line):
                break;
            row = line.split("=")
            key = row[0].strip()
            val = row[1].strip()
            self.namelist[key] = val
            self[key] = val

    def __setattr__(self, name, value):
        super(Main, self).__setattr__(name, value)

    def __setitem__(self, name, value):
        self.__setattr__(name, value)

    def __getattribute__(self, name):
        return(object.__getattribute__(self, name))

    def __getitem__(self, name):
        return(self.__getattribute__(name))

def qm_handler(func):
    def wrapper_qm_choice(**kwargs):
        return(func(**kwargs))
    return(wrapper_qm_choice)

@qm_handler
def qm_choice(infile = None, tpl = None):
    """Return a QM object, depending on the template used"""
    if(tpl == None):
        return(QMtemplate())
    elif(tpl == "nwchem"):
        return(NWChemtpl(infile))
    elif(tpl == "gaussian"):
        return(Gautpl(infile))
    elif(tpl == "molcas"):
        return(Molcastpl(infile))
    else:
        return("No function returned")

def geometry_handler(traj, mainobj, tplobj):
    """This function generates one input file per geometry requested.
    If a range of geometries was chosen it generates the inputs sequentially. 
    Consider also calling the main script in parallel.
    Args: 
        traj    = pytraj trajectory file
        mainobj = a Main() object
        tplobj  = a QMtemplate() object 
    """
    currdir = os.getcwd() + "/"
    igeom = int(mainobj.geoms)
    _geometry_handler(igeom, traj, mainobj, tplobj, currdir)
#    has_geoms = (mainobj.geoms != None)
#    if(has_geoms):
#        # Iterate over the range of geometries
#        geoms = mainobj.geoms
#        rng = geoms.split()
#        if(len(rng) == 3):
#            start, stop, step = [int(igeom) for igeom in rng[:3]] 
#        elif(len(rng) == 2):
#            start, stop = [int(igeom) for igeom in rng[:2]] 
#            step = 1
#        elif(len(rng) == 1):
#            start = int(rng[0])
#            stop  = start
#            step  = 1
#        else:
#            raise IOError("Select a valid number of geometries in geoms")
#        stop += 1
#
#        for igeom in range(start, stop, step):
#            # Generate directory for ith geometry
#            _geometry_handler(igeom, traj, mainobj, tplobj, currdir)
#            os.chdir(currdir)
#    else:
#        igeom = int(mainobj.igeom)
#        _geometry_handler(igeom, traj, mainobj, tplobj, currdir)
#        os.chdir(currdir)
    os.chdir(currdir)

def _geometry_handler(igeom, traj, mainobj, tplobj, currdir):
    """Executed from within geometry_handler"""
    idir = "{:s}geom{:d}/".format(currdir, igeom)
    try:
        os.mkdir(idir)
    except:
        pass
    os.chdir(idir)
    _closest_handler(traj[igeom:igeom+1], mainobj, tplobj)
    tplobj.write_input(traj, traj.top, mainobj.qmmask, igeom = igeom)

# Assess whether to analyze closest residues or not
def _closest_handler(traj, mainobj, tplobj):
    """Handle situation in which a number of closest
    solvent molecules is required"""
    has_solvent = mainobj.solvmask != None
    has_closest = mainobj.closest != None
    has_bsse    = tplobj.bsse != None
    if(has_solvent & has_closest):
        print("Analyzing the {:s} closest residues to the QM mask {}\n".format(mainobj.closest, mainobj.qmmask) )
        clmasks = get_closest(traj, mainobj)
        # Update main.qmmask attribute, with the extended qmmask
        mainobj["qmmask"] = clmasks[0]
        # Define bsse attribute for tplobj, if requested
        if(has_bsse):
            tplobj["bsse"] = ["mon{:d} = {:s}".format(i+1, imask) for i, imask in enumerate(clmasks[1:])]
#            print("tplobj.bsse = ", tplobj.bsse)

    else:
        print("No closest feature requested")
        pass


def get_closest(traj, mainobj):
    """Return qmmask for the complex, the ligand (mon1) and the environment\
            (mon2) considering the selected number of closest residues"""
    qmmask   = mainobj.qmmask
    solvent  = mainobj.solvmask
    closest  = mainobj.closest

    # Get lists
    residues = _get_complex_list(traj, qmmask, solvent, closest)
    ligand   = _get_ligand_list(traj, qmmask)
    env      = list(np.setdiff1d(residues, ligand))
    # Get masks
    cplmask  = _get_mask_res(residues) 
    ligmask  = _get_mask_res(ligand) 
    envmask  = _get_mask_res(env)
    return(cplmask, ligmask, envmask)

def _get_mask_res(mask_list):
    """Get a compact mask in terms of residues"""
    mlist = list(map(str, mask_list))
    mask = ":" + ",".join(mlist)
    return(mask)

def _get_ligand_list(traj, qmmask):
    """To be used within get_closest. This function returns a list
    of the residues present in qmmask (the "ligand" mask).
    """
    ligres = []
    qmtraj = traj[qmmask]
    for ires in traj.top.residues:
        for iatom in range(ires.first_atom_index, ires.last_atom_index):
            icrd = list(traj.xyz[0][iatom])
            if(icrd in qmtraj.xyz[0]):
                atom = traj.top.atom(iatom)
                ligres.append(atom.resid + 1) # Indexing starts from zero
    ligres = list(np.unique(ligres))
    return(ligres)

def _get_complex_list(traj, qmmask, solvent, closest):
    """To be used within get_closest. This function returns a list
    of the closest N residues present around qmmask (the "complex" mask).
    """
    residues = []
    mask = "{},{}".format(qmmask, solvent)
#    print("mask in _get_complex_list = ", mask)
    cltrj    = pt.closest(traj[mask], mask = qmmask, solvent_mask = solvent,\
            n_solvents = closest, dtype = "trajectory")
    
    # Iterate over all residues in traj, then retrieve all indices of those
    # residues belonging to cltrj
    for ires in traj.top.residues:
        for iatom in range(ires.first_atom_index, ires.last_atom_index):
            icrd = traj.xyz[0][iatom]
            # Element-wise comparison
            assess = [list(icrd) == list(cltrj.xyz[0][i]) for i in range(len(cltrj.xyz[0]))]
            if(True in assess):
                atom = traj.top.atom(iatom)
                residues.append(atom.resid + 1) # Indexing starts from zero
    residues = list(np.unique(residues))
    print("residues in the QM region = {}\n".format(residues))
    return(residues)


def range_geometries(init_main, mainfile, templatefile):
    """Iterate over a range of geometries"""
    has_geoms = (init_main.geoms != None)

    if(has_geoms):
        # Iterate over the range of geometries
        geoms = init_main.geoms
        # Define template object
        rng = geoms.split()
        if(len(rng) == 3):
            start, stop, step = [int(igeom) for igeom in rng[:3]] 
        elif(len(rng) == 2):
            start, stop = [int(igeom) for igeom in rng[:2]] 
            step = 1
        elif(len(rng) == 1):
            start = int(rng[0])
            stop  = start
            step  = 1
        else:
            raise IOError("Select a valid number of geometries in geoms")
        stop += 1

        for igeom in range(start, stop, step):
            # Generate directory for ith geometry
            # Define mainobj for ith geometry
            print("#################### GENERATING INPUT FILE FOR GEOMETRY {} ####################\n".format(igeom))
            mainobj = Main(mainfile)
            mainobj.geoms = igeom
            traj = pt.load(mainobj.traj, mainobj.top)
            tplobj = qm_choice(tpl = mainobj.tpl, infile = templatefile)
            geometry_handler(traj, mainobj, tplobj)
            print("#################### END OF INPUT GENERATION               ####################\n\n".format(igeom))
    else:
        raise IOError("Select a valid number of geometries in geoms")
#        igeom = int(mainobj.igeom)
#        geometry_handler(traj, mainobj, tplobj)
#        os.chdir(currdir)



if(__name__ == "__main__"):

    parser = ArgumentParser("""MoBioTools QM/MM input generator""")
    parser.add_argument("-i", dest = "mainfile", type = str, help = "Input file")
    parser.add_argument("-t", dest = "tplfile", type = str, help = "Template file")

    options   = parser.parse_args()
    mainfile  = options.mainfile
    template  = options.tplfile
    
    init_main = Main(mainfile)

    range_geometries(init_main, mainfile, template)
#    qmtpl = qm_choice(tpl = main.tpl, infile = template)

    # Load trajectory
#    traj = pt.load(main.traj, main.top)
#    geometry_handler(traj, main, qmtpl)
    print("If you are seeing this, it means the program has attained a Normal Termination. A Tope!")
    
