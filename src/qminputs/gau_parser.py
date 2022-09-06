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

class Gautpl(QMtemplate):

    def __init__(self, infile):
        """Parse Gaussian template file"""
        QMtemplate.__init__(self)

        # Attrubutes specific to gaussian
        self._infile     = infile
        self._template   = "gaussian"
        self.namelist.update({"route": None, "scrf": None})

        # Parse common sections
        self.parse(infile)
    
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
                    igeom = 0, 
                    charge = 0, 
                    spin = 1,
                    null_charges = False,
                    prename = None,
                    bqmask = None):
        """traj = pytraj trajectory, top = pytraj topology"""
        self.traj    = traj
        self.top     = top
        self.qmmask  = qmmask
        self.chgmask =  "!(" + qmmask + ")"

        # Define variables and defaults
        if (prename != None):
            outfile = prename + "_geom" + str(igeom) + ".com"
        else:
            outfile = "geom" + str(igeom) + ".com"
        
        # Define link handler
        self.link_handler(traj, qmmask, igeom)

        # Write output
        with open(outfile, "w") as f:
            igeom = int(igeom)
            self._write_header(f, igeom, outfile)
            self._write_route(f)
            self._write_title(f, outfile)
            self._write_chgspin(f, charge, spin)
            self._write_geometry(f, igeom, bqmask = bqmask)
            self._write_point_charges(f, igeom, null_charges)
            self._write_basis(f)
            self._write_final_space(f)
    
            for attr in self.namelist.keys():
                self._write_other(f, attr)

        # Write xyz for QM and MM regions (for 
        # visualization purposes)
        self.linkobj.write_geometry(True)

    def _write_final_space(self, f):
        """"""
        f.write("\n")

    def _write_header(self, f, igeom, outfile = None):
        """Write Link0 commands"""
        for iopt in self.namelist["header"]:
            if("chk" in iopt.lower()):
                f.write("%chk=" + outfile.replace(".com", ".chk") + "\n")
            else:
                f.write(iopt + "\n")

    def _write_route(self, f):
        """Write route section"""
        try:
            f.write(self.namelist["route"][0] + "\n")
        except:
            print("No route section found")

    def _write_title(self, f, title):
        f.write("\n" + title + "\n")

    def _write_chgspin(self, f, charge = 0, spin = 1):
        """Write charge and spin headers"""
        has_frag = self.namelist["fragments"] != None
        if(has_frag):
            f.write("\n" + self.chgspin[0] + "\n")
        else:
            f.write("\n" + str(charge) + " " + str(spin) + "\n")

    def _write_geometry(self, f, igeom, bqmask = None):
        """Write xyz coordinates. f = file opened with the append option"""

        fmt = "{:<7s}{:<15.6f}{:<15.6f}{:<15.6f}\n"
        # Iterator containing all atoms in the QM region
        itatoms = self.traj[self.qmmask].top.atoms

        if(bqmask != None):
            bqcrd = []
            # Get bqmask atom coordinates
            for ibqc in self.traj[bqmask].xyz[igeom]:
                bqcrd.append(list(ibqc))

            for cnt, icrd in enumerate(self.traj[self.qmmask].xyz[igeom]):
                if(list(icrd) in bqcrd):
                    iatom = ATOMNUM[next(itatoms).atomic_number].upper() + "-Bq"
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
                    link_at = self.linkobj.link_atom + "-Bq"
                else:
                    link_at = self.linkobj.link_atom
                f.write(fmt.format(link_at, *icrd))

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
        f.write("\n")

    def _write_point_charges(self, f, igeom, null_charges = False):
        """Write external point charges. 
        null_charges = whether to set to zero the point charges 
        (e.g: in the case of a ligand+environment system)Consider eventual link atoms.
        """
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
                for nat, iat in enumerate(self.traj.top.atoms):
                    if(self.linkobj.labels[nat] == 1):
                        # Consider only MM atoms
                        if(nat not in self.linkobj.all_link):
                            # MM atoms not involved in QM/MM boundary
                            crd = self.traj.xyz[igeom][nat]
                            chg = charges[nat]
                            f.write(fmt.format(*crd, chg))
                        elif(nat in self.linkobj.qmmm_bonds[:,1]):
                            # Print NNs for each QM/MM atom with smeared charges
                            chg_parz = charges[nat]/\
                                       len(self.linkobj.nn_atoms[nat])
                            
                            for inn in self.linkobj.nn_atoms[nat]:
                                # Ith nearest neighbor
                                crd = self.traj.xyz[igeom][inn]
                                chg = charges[inn] + chg_parz
                                f.write(fmt.format(*crd, chg))
                        else:
                            pass
            f.write("\n")
        else:
            pass
        
    def _write_other(self, f, attribute):
        """Write other nwchem sections"""
#        exclude  = ["header", "chgspin", "route", "externchg", "bsse"]
        exclude  = ["header", "chgspin", "route", "externchg", "bsse", "basis"]
        keys     = np.setdiff1d(list(self.namelist.keys()), exclude)
#        print("other keys = ", keys)
        has_attr = (attribute in keys and self.namelist[attribute] != None)

        if(has_attr):
#            f.write("\n")
            for iopt in self.namelist[attribute]:
                f.write(iopt + "\n")
            f.write("\n")
    
    def _write_basis(self, f):
        """Write basis at the bottom if the /gen keyword is evidenced on the route"""
        if("/gen" in self.route[0].lower()):
            for iopt in self.namelist["basis"]:
                f.write(iopt + "\n")
            f.write("\n")

