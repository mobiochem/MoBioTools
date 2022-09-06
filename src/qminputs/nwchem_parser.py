#!/usr/bin/env python3
import numpy as np
import pytraj as pt
from copy import deepcopy
from collections import OrderedDict
from argparse import ArgumentParser
import os
from constants import ATOMNUM
from common_parser import QMtemplate, write_input_bsse 
from link_atoms import Link_atoms, get_unit_vector

class NWChemtpl(QMtemplate):

    def __init__(self, infile):
        """Parse NWChem template file"""
        QMtemplate.__init__(self)
        
        # Attrubutes specific to nwchem
        self._infile     = infile
        self._template   = "nwchem"
        self.namelist.update({"dft": None, "tasks": None, "cosmo": None})

        # Parse common sections
        self.parse(self._infile)

    def link_handler(self, traj, qmmask, igeom = 0):
        """Initialize Link_atom object to handle
           eventual link atoms
        """

        self.linkobj = Link_atoms(traj, qmmask, igeom)
        self.linkobj.get_link_info()
        self.N_link = self.linkobj.N_link
        if(self.N_link>0):
            print("QM/MM boundaries found: printing bond matrix\n\
                   involving those bonds in the format \n\
                   [[qm_id1, mm_id1],\n\
                    [qm_id2, mm_id2], ...]")
            print(self.linkobj.qmmm_bonds)
            print("Employing link atom approach...")
        else:
            print("No QM/MM boundary found")

    @write_input_bsse
    def write_input(self, 
                    traj, 
                    top, 
                    qmmask, 
                    igeom        = 0, 
                    charge       = 0, 
                    spin         = 1,
                    null_charges = False,
                    prename      = None,
                    bqmask       = None):
        """traj = pytraj trajectory, top = pytraj topology"""
        self.traj    = traj
        self.top     = top
        self.qmmask  = qmmask
        self.chgmask =  "!(" + qmmask + ")"
        
        # Define variables and defaults
        if (prename != None):
            outfile = prename + "_geom" + str(igeom) + ".in"
        else:
            outfile = "geom" + str(igeom) + ".in"

        # Set spin on dft attribute
        has_dft = self.namelist["dft"] != None
        if(has_dft):
            id_mult = np.where(["mult" in iop for iop in self.namelist["dft"]])[0][0]
            print("dft namelist", self.namelist["dft"])
            print("id mult = ", id_mult)
            self.namelist["dft"][id_mult] = "mult " + str(spin)
        
        # Define link handler
        self.link_handler(traj, qmmask, igeom)
        
        # Write output
        with open(outfile, "w") as f:
            igeom = int(igeom)
            self._write_header(f, igeom)
            self._write_geometry(f, igeom, bqmask = bqmask)
            self._write_charge(f, charge)
            self._write_point_charges(f, igeom, prename, null_charges)
            
            for attr in self.namelist.keys():
                self._write_other(f, attr)
            self._write_tasks(f)
        
        # Write xyz for QM and MM regions (for 
        # visualization purposes)
        self.linkobj.write_geometry(True)

    def _write_header(self, f, igeom):
        for iopt in self.namelist["header"]:
            if(iopt.split()[0].lower() == "start"):
                f.write("start geom" + str(igeom) + "\n")
            elif(iopt.split()[0].lower() == "title"):
                f.write("title geom" + str(igeom) + ".in\n")
            else:
                f.write(iopt + "\n")
#        f.write("end\n")

    def _write_geometry(self, f, igeom, symmetry = None, center = None, bqmask = None):
        """Write xyz coordinates. f = file opened with the append option"""
        
        fmt = "{:<7s}{:<15.6f}{:<15.6f}{:<15.6f}\n"
        # Iterator containing all atoms in the QM region
        itatoms = self.traj[self.qmmask].top.atoms
        f.write("geometry noautosym nocenter\n") # Avoid displacements. Also no symmetry

        # BQMask
        if(bqmask != None):
            bqcrd = []
            # Get bqmask atom coordinates
            for ibqc in self.traj[bqmask].xyz[igeom]:
                bqcrd.append(list(ibqc))

            for cnt, icrd in enumerate(self.traj[self.qmmask].xyz[igeom]):
                if(list(icrd) in bqcrd):
                    iatom = "bq" + ATOMNUM[next(itatoms).atomic_number].upper()
                    f.write(fmt.format(iatom, *icrd))
                else:
                    iatom = ATOMNUM[next(itatoms).atomic_number].upper()
                    f.write(fmt.format(iatom, *icrd))
            # Now print eventual Link atoms, accounting for
            # those present in the Bq region
            for cnt, icrd in enumerate(self.linkobj.link_pos):
                qm_id  = self.linkobj.qmmm_bonds[cnt][0]
                crd_qm = self.traj.xyz[igeom][qm_id] 
                # Assess whether the QM partner is in the Bq region
                if(list(crd_qm) in bqcrd):
                    link_at = "bq" + self.linkobj.link_atom
                else:
                    link_at = self.linkobj.link_atom
                f.write(fmt.format(link_at, *icrd))
            f.write("end\n")
        
        # No Bq atoms present
        else:
            for cnt, icrd in enumerate(self.traj[self.qmmask].xyz[igeom]):
                iatom = ATOMNUM[next(itatoms).atomic_number].upper()
                f.write(fmt.format(iatom, *icrd))
            # Now print eventual Link atoms
            for cnt, icrd in enumerate(self.linkobj.link_pos):
                qm_id   = self.linkobj.qmmm_bonds[cnt][0]
                crd_qm  = self.traj.xyz[igeom][qm_id] 
                link_at = self.linkobj.link_atom
                f.write(fmt.format(link_at, *icrd))
            f.write("end\n")

    def _write_charge(self, f, charge = 0):
        """Write global charge"""
        f.write("charge " + str(charge) + "\n")

    def _write_point_charges(self, f, igeom, prename = None, null_charges = False):
        """Write external point charges
        print_chg = whether to generate a point charges file
        null_charges     = if True, set to zero the point charges"""
#        try:
#            print_chg = (self.namelist["externchg"][0] == "true")
#        except:
#            print("No point charges evidenced. Do not generate a point charge file")
#        print_chg = False
        print_chg = self.namelist["externchg"] != None
        print("print_chg = ", print_chg)
        if(print_chg):
            has_chg = len(self.traj[self.chgmask])>0
            if(has_chg):
                if(null_charges):
                    charges = np.zeros(len(self.traj.top.charge))
                else:
                    charges = self.traj.top.charge
                fmt = "{:<15.6f}{:<15.6f}{:<15.6f}{:<15.6f}\n"
                
                # Write bq option on input file (for point charges)
                if(prename != None):
                    chgfile = "charges_{:s}_geom{:d}.xyz".format(prename,igeom)
                else:
                    chgfile = "charges_geom{:d}.xyz".format(igeom)
                f.write("bq\n")
                f.write("load " + chgfile + " format 1 2 3 4\n")
                f.write("end\n")
                
                # Write charges file
                with open(chgfile, "w") as g:
                    for nat, iat in enumerate(self.traj.top.atoms):
                        if(self.linkobj.labels[nat] == 1):
                            # Consider only MM atoms
                            if(nat not in self.linkobj.all_link):
                                # MM atoms not involved in QM/MM boundary
                                crd = self.traj.xyz[igeom][nat]
                                chg = charges[nat]
                                g.write(fmt.format(*crd, chg))
                            elif(nat in self.linkobj.qmmm_bonds[:,1]):
                                # Print NNs for each QM/MM atom with smeared charges
                                chg_parz = charges[nat]/\
                                           len(self.linkobj.nn_atoms[nat])
                                
                                for inn in self.linkobj.nn_atoms[nat]:
                                    # Ith nearest neighbor
                                    crd = self.traj.xyz[igeom][inn]
                                    chg = charges[inn] + chg_parz
                                    g.write(fmt.format(*crd, chg))
                            else:
                                pass

        else:
            pass

    def _write_other(self, f, attribute):
        """Write other nwchem sections"""
        exclude  = ["header", "chgspin", "tasks", "externchg", "bsse"]
        keys     = np.setdiff1d(list(self.namelist.keys()), exclude)
        print("other keys = ", keys)
        has_attr = (attribute in keys and self.namelist[attribute] != None)

        if(has_attr):
            f.write(attribute + "\n")
            for iopt in self.namelist[attribute]:
                f.write(iopt + "\n")
            f.write("end\n")

    def _write_tasks(self, f):
        """Write tasks"""
        has_attr = "tasks" in self.namelist.keys()
        if(has_attr):
            for iopt in self.namelist["tasks"]:
                f.write(iopt + "\n")

