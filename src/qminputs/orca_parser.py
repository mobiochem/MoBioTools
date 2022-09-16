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

    def link_handler(self, traj, qmmask, igeom = 0, link_atom = "H", link_dist = 1.09):
        """Initialize Link_atom object to handle
           eventual link atoms
        """

        self.linkobj = Link_atoms(traj, qmmask, igeom, link_atom, link_dist)
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
                    bqmask = None,
                    link_atom = "H",
                    link_dist = 1.09):
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
        
        # Define link handler
        self.link_handler(traj, qmmask, igeom, link_atom, link_dist)

        # Write output
        with open(outfile, "w") as f:
            igeom = int(igeom)
            self._write_header(f)
            self._write_geometry(f, igeom, charge, spin, bqmask = bqmask)
            self._write_point_charges(f, igeom, prename, null_charges)

            for attr in self.namelist.keys():
                self._write_other(f, attr)

        # Write xyz for QM and MM regions (for 
        # visualization purposes)
        self.linkobj.write_geometry(True)

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

        fmt = "{:<7s}{:<15.6f}{:<15.6f}{:<15.6f}\n"
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
                    link_at = self.linkobj.link_atom + ":"
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
                    charges = np.zeros(len(self.traj.top.charge))
                else:
                    charges = self.traj.top.charge
                fmt = "{:<15.6f}{:<15.6f}{:<15.6f}{:<15.6f}\n"
                
                # Write %pointcharges option on input file (for point charges)
                if(prename != None):
                    chgfile = "charges_{:s}_geom{:d}.xyz".format(prename,igeom)
                else:
                    chgfile = "charges_geom{:d}.xyz".format(igeom)
                f.write('%pointcharges "{}"\n'.format(chgfile))
                
                # Write charges file
                with open(chgfile, "w") as g:
                    MM_charges = self.linkobj.N_mm - self.linkobj.N_link

                    g.write("{:d}\n\n".format(MM_charges))
#                    g.write("{:d}\n\n".format(len(charges)))
#                    for cnt, icrd in enumerate(self.traj[self.chgmask].xyz[igeom]):
#                        ichg = charges[cnt]
#                        g.write(fmt.format(ichg, *icrd) + "\n")
                    for nat, iat in enumerate(self.traj.top.atoms):
                        if(self.linkobj.labels[nat] == 1):
                            # Consider only MM atoms
                            if(nat not in self.linkobj.all_link):
                                # MM atoms not involved in QM/MM boundary
                                crd = self.traj.xyz[igeom][nat]
                                chg = charges[nat]
                                g.write(fmt.format(chg, *crd))
                            elif(nat in self.linkobj.qmmm_bonds[:,1]):
                                # Print NNs for each QM/MM atom with smeared charges
                                chg_parz = charges[nat]/\
                                           len(self.linkobj.nn_atoms[nat])
                                
                                for inn in self.linkobj.nn_atoms[nat]:
                                    # Ith nearest neighbor
                                    crd = self.traj.xyz[igeom][inn]
                                    chg = charges[inn] + chg_parz
                                    g.write(fmt.format(chg, *crd))
                            else:
                                pass
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



