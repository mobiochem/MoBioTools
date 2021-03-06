#!/usr/bin/env python3
import numpy as np
import pytraj as pt
from copy import deepcopy
from collections import OrderedDict
from argparse import ArgumentParser
import os
from constants import ATOMNUM
from common_parser import QMtemplate, write_input_bsse 

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

        fmt = "{:<7s}{:<15.6f}{:<15.6f}{:<15.6f}"
        # Iterator containing all atoms in the topology file
#        itatoms = self.top.atoms
        itatoms = self.traj[self.qmmask].top.atoms

        if(bqmask != None):
            bqcrd = []
            # Get bqmask atom coordinates
            for ibqc in self.traj[bqmask].xyz[igeom]:
                bqcrd.append(list(ibqc))

            for cnt, icrd in enumerate(self.traj[self.qmmask].xyz[igeom]):
                if(list(icrd) in bqcrd):
                    iatom = ATOMNUM[next(itatoms).atomic_number].upper() + "-Bq"
                    f.write(fmt.format(iatom, *icrd) + "\n")
                else:
                    iatom = ATOMNUM[next(itatoms).atomic_number].upper()
                    f.write(fmt.format(iatom, *icrd) + "\n")
        else:
            for cnt, icrd in enumerate(self.traj[self.qmmask].xyz[igeom]):
                iatom = ATOMNUM[next(itatoms).atomic_number].upper()
                f.write(fmt.format(iatom, *icrd) + "\n")
        f.write("\n")


    def _write_point_charges(self, f, igeom, null_charges = False):
        """Write external point charges. 
        null_charges = whether to set to zero the point charges 
        (e.g: in the case of a ligand+environment system)"""
#        try:
#            print_chg = (self.namelist["externchg"][0] == "true")
#        except:
#            print("No point charges evidenced. Do not generate point charges")
#            print_chg = False
        print_chg = self.namelist["externchg"] != None
        print("print_chg = ", print_chg)
        if(print_chg):
            has_chg = len(self.traj[self.chgmask])>0
            if(has_chg):
                if(null_charges):
                    charges = np.zeros(len(self.traj[self.chgmask].top.charge))
                else:
                    charges = self.traj[self.chgmask].top.charge
                fmt = "{:<15.6f}{:<15.6f}{:<15.6f}{:<15.6f}"
                # Write point charges
                for cnt, icrd in enumerate(self.traj[self.chgmask].xyz[igeom]):
                    ichg = charges[cnt]
                    f.write(fmt.format(*icrd, ichg) + "\n")
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

