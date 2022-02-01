#!/usr/bin/env python3
import numpy as np
import pytraj as pt
from copy import deepcopy
from collections import OrderedDict
from argparse import ArgumentParser
import os
from constants import ATOMNUM
from common_parser import QMtemplate, write_input_bsse 

class Orcatpl(QMtemplate):

    def __init__(self, infile):
        """Parse Orca template file"""
        QMtemplate.__init__(self)

        # Attrubutes specific to gaussian
        self._infile     = infile
        self._template   = "orca"
        self.namelist.update({"tddft": None, 
                              "xtb": None, 
                              "casscf": None, 
                              "cpcm": None,
                              "casscf": None})

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
            outfile = prename + "_geom" + str(igeom) + ".inp"
        else:
            outfile = "geom" + str(igeom) + ".inp"

        # Write output
        with open(outfile, "w") as f:
            igeom = int(igeom)
            self._write_header(f)
            self._write_geometry(f, igeom, charge, spin, bqmask = bqmask)
            self._write_point_charges(f, igeom, prename, null_charges)

            for attr in self.namelist.keys():
                self._write_other(f, attr)

    def _write_header(self, f):
        """Write KEYWORD commands (introduced by !)"""
        for iopt in self.namelist["header"]:
            if("!" not in iopt):
                f.write("! {}\n".format(iopt))
            else:
                f.write("{}\n".format(iopt))
        f.write("\n")

    def _write_geometry(self, f, igeom, charge = 0, spin = 1, bqmask = None):
        """Write xyz coordinates. They are written in conjunction with the
        charge and the spin multiplicity. The output format is the following:
        * xyz charge spin
        H x_coord y_coord z_coord
        ...
        *
        """

        fmt = "{:<7s}{:<15.6f}{:<15.6f}{:<15.6f}"
        itatoms = self.traj[self.qmmask].top.atoms
        # Write section header
        f.write("* xyz {} {}\n".format(charge, spin))
        if(bqmask != None):
            bqcrd = []
            # Get bqmask atom coordinates
            for ibqc in self.traj[bqmask].xyz[igeom]:
                bqcrd.append(list(ibqc))

            for cnt, icrd in enumerate(self.traj[self.qmmask].xyz[igeom]):
                if(list(icrd) in bqcrd):
                    iatom = ATOMNUM[next(itatoms).atomic_number].upper() + " : "
                    f.write(fmt.format(iatom, *icrd) + "\n")
                else:
                    iatom = ATOMNUM[next(itatoms).atomic_number].upper()
                    f.write(fmt.format(iatom, *icrd) + "\n")
        else:
            for cnt, icrd in enumerate(self.traj[self.qmmask].xyz[igeom]):
                iatom = ATOMNUM[next(itatoms).atomic_number].upper()
                f.write(fmt.format(iatom, *icrd) + "\n")
        # Write section closure
        f.write("*\n\n")

    def _write_point_charges(self, f, igeom, prename = None, null_charges = False):
        """Write external point charges. 
        null_charges = whether to set to zero the point charges 
        (e.g: in the case of a ligand+environment system)"""
        
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
                
                # Write %pointcharges option on input file (for point charges)
                if(prename != None):
                    chgfile = "charges_{:s}_geom{:d}.xyz".format(prename,igeom)
                else:
                    chgfile = "charges_geom{:d}.xyz".format(igeom)
                f.write('%pointcharges "{}"\n'.format(chgfile))
                
                # Write charges file
                with open(chgfile, "w") as g:
                    g.write("{:d}\n\n".format(len(charges)))
                    for cnt, icrd in enumerate(self.traj[self.chgmask].xyz[igeom]):
                        ichg = charges[cnt]
                        g.write(fmt.format(*icrd, ichg) + "\n")
                f.write("\n")
        else:
            pass

    def _write_other(self, f, attribute):
        """Write other orca sections introduced by % and 
        closed by %end
        """
        exclude  = ["header", "chgspin", "tasks", "externchg", "bsse"]
        keys     = np.setdiff1d(list(self.namelist.keys()), exclude)
        print("other keys = ", keys)
        has_attr = (attribute in keys and self.namelist[attribute] != None)

        if(has_attr):
            f.write("%{}\n".format(attribute))
            for iopt in self.namelist[attribute]:
                f.write(iopt + "\n")
            f.write("end\n")



