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

class Molcastpl(QMtemplate):

    def __init__(self, infile):
        """Parse Molcas template file"""
        QMtemplate.__init__(self)

        # Attrubutes specific to gaussian
        self._infile     = infile
        self._template   = "molcas"
        self.namelist.update({"seward": None, 
                              "scf": None,
                              "rasscf": None,
                              "caspt2": None,
                              "rassi": None,
                              "guess": None})

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

        print("qmmask for prename", prename, " = ", qmmask)
        # Define variables and defaults
        if (prename != None):
            outfile = prename + "_geom" + str(igeom) + ".in"
        else:
            outfile = "geom" + str(igeom) + ".in"

        # Define link handler
        self.link_handler(traj, qmmask, igeom, link_atom, link_dist)

        # Prepare guess files, if requested
        # No guess for SCF
        guess_tpl = {"inporb": outfile + ".guess.RasOrb",
                     "jobiph": outfile + ".guess.JobIph",
                     "jobmix": outfile + ".guess.JobMix"}
        guess_dict = {"inporb": None,
                      "jobiph": None,
                      "jobmix": None}

        if(self.namelist["guess"] != None):
            guess_dict = {key: guess_tpl[key] for key in guess_tpl.keys() \
                    if key in self.guess}
        
        with open(outfile, "w") as f:
            igeom = int(igeom)
            self._write_header(f, igeom, prename)
            self._write_seward(f)
            self._write_geometry(igeom, bqmask, prename)
            self._write_point_charges(f, igeom, prename, null_charges)
            # Write (eventually) SCF, RASSCF, CASPT2 and RASSI
            if(self.namelist["scf"] != None ):
                # No guess for SCF
                self._write_calculation(f, calctype = "scf", 
                                          outfile = outfile)
#                                          inporb  = inporb)
#                inporb = outfile."ScfOrb"
                guess_dict["inporb"] = outfile + ".ScfOrb"
            if(self.rasscf != None):
                self._write_calculation(f, calctype = "rasscf", 
                                          outfile = outfile, 
                                          inporb  = guess_dict["inporb"])
#                                          inporb  = inporb)
#                inporb = outfile."RasOrb"
#                jobiph = outfile."JobIph"
                guess_dict["inporb"] = outfile + ".RasOrb"
                guess_dict["jobiph"] = outfile + ".JobIph"
            if(self.caspt2 != None):
                self._write_calculation(f, calctype = "caspt2", 
                                          outfile = outfile, 
                                          jobiph  = guess_dict["jobiph"])
#                jobmix = outfile."JobMix"
                guess_dict["jobmix"] = outfile + ".JobMix"
            if(self.rassi != None):
                self._write_calculation(f, calctype = "rassi", 
                                          outfile = outfile, 
                                          jobmix  = guess_dict["jobmix"])
#                jobmix = outfile."JobMix"
        
        # Write xyz for QM and MM regions (for 
        # visualization purposes)
        self.linkobj.write_geometry(True)
          
    def _write_seward(self, f):
        """Write seward command"""
        if(self.seward != None):
            f.write("&SEWARD &END\n")

    def _write_header(self, f, igeom, prename = None):
        """Write gateway commands"""

        # Copy point charges, if required
        """Add point charges file name in header section"""
        """Define chgfile if point charges are to be copied"""
        if(self.namelist["externchg"] != None):
            if(self.bsse != None and prename != None):
                chgfile = "charges_{:s}_geom{:d}.xyz".format(prename, igeom)
            else:
                chgfile = "charges_geom{:d}.xyz".format(igeom)
            self._copy_to_workdir(f, chgfile, chgfile)
        else:
            chgfile = None


        f.write("&GATEWAY &END\n")
        
        """Add name of xyz file in gateway (Two if bsse)"""
        if(self.bsse != None and prename != None and prename != "complex"):
            mons  = ["mon1", "mon2"]
#            for imon in self.bsse:
            for imon in mons:
                namefile = "coord = {:s}_geom{:d}.xyz\n".format(imon, igeom)
                f.write(namefile)
#                self.header.append(namefile)
            # Print the assigned bq mask
            index = np.where(np.array(mons) == prename)[0]
            bqidx = np.setdiff1d([0,1], index)[0]
            f.write("bsse = {:d}\n".format(bqidx + 1))
        else:
            namefile = "coord = geom{:d}.xyz\n".format(igeom)
            f.write(namefile)
#            self.header.append(namefile)
        
        """Add point charges file name in header section, if required"""
        if(chgfile != None):
            f.write("Xfield = {:s}\n".format(chgfile))
#        if(self.namelist["externchg"] != None):
#            if(self.bsse != None and prename != None):
#                chgfile = "charges_{:s}_geom{:d}.xyz".format(prename, igeom)
#            else:
#                chgfile = "charges_geom{:d}.xyz".format(igeom)
#            f.write("Xfield = {:s}\n".format(chgfile))
#            # 24/03/2022: Copy point charges file from CurrDir to WorkDir
#            self._copy_to_workdir(f, chgfile, chgfile)

        for iopt in self.namelist["header"]:
            f.write(iopt + "\n")

        
    
    def _write_calculation(self, f, calctype = None, 
                            outfile  = None, 
                            inporb   = None, 
                            jobiph   = None,
                            jobmix   = None):
        """Write Calculation section;
           calctype = SCF, RASSCF, CASPT2, RASSI
           outfile  = MolCas input file name
           inporb   = input orbitals file name
           jobiph   = guess wavefunction file name
           jobmix   = guess for rassi
        """

        if(self[calctype.lower()] != None):
            # Assess whether there are guess files to copy
            guess = {inporb: "INPORB", jobiph: "JOBIPH", jobmix: "JOBMIX"}
            for iguess in guess.keys():
                if(iguess != None):
                    self._copy_to_workdir(f, iguess, guess[iguess])
            # Write header section
            f.write("&{:s} &END\n".format(calctype.upper()))
#            for iopt in self[calctype.lower()]:
#                f.write(iopt + "\n")
            molden = outfile + ".{:s}.molden".format(calctype.lower())
            mos    = outfile + ".{:s}Orb".format(calctype[0].upper() + calctype[1:3].lower())
            wfn    = outfile + ".JobIph"
            wfnmix = outfile + ".JobMix"
            
            # Write contents, if any
            list_calc = self[calctype.lower()]
            if(len(list_calc)>0):
                for irow in list_calc:
                    f.write(irow + "\n")
            
            if(calctype.lower() == "scf" or calctype.lower() == "rasscf") :
                self._copy_to_currdir(f, molden ,molden)
                self._copy_to_currdir(f, mos ,mos)
#            elif(calctype.lower() == "caspt2"):
#                self._copy_to_currdir(f, wfn, wfn)
            elif(calctype.lower() == "caspt2"):
                self._copy_to_currdir(f, wfnmix, wfnmix)
            elif(calctype.lower() == "rassi"):
                self._copy_to_currdir(f, wfnmix, wfnmix)



    def _copy_to_currdir(self, f, source_file, target_file):
        """Print copy command on the input file"""
        f.write(">>COPY $WorkDir/{:s} $CurrDir/{:s}\n".format(source_file, target_file))

    def _copy_to_workdir(self, f, source_file, target_file):
        """Print copy command on the input file"""
        f.write(">>COPY $CurrDir/{:s} $WorkDir/{:s}\n".format(source_file, target_file))

    def _write_geometry(self, igeom, bqmask = None, prename = None):
        """Write xyz coordinates. If bqmask != None, generate two xyz
        files, each one containing one of the monomers. In this case the 
        function prints out the atoms NOT present in bqmask."""

        fmt = "{:<7s}{:<15.6f}{:<15.6f}{:<15.6f}\n"
        # Iterator containing all atoms in the topology file
#        itatoms = self.top.atoms
        itatoms = self.traj[self.qmmask].top.atoms

        if(bqmask != None and prename != None):
            bqcrd = []
            outfile = "{:s}_geom{:d}.xyz".format(prename, igeom)
            with open(outfile, "w") as f:
                # Get bqmask atom coordinates
                for ibqc in self.traj[bqmask].xyz[igeom]:
                    bqcrd.append(list(ibqc))

                # Print number of atoms that are not in bqmask
                mask = "(" + self.qmmask + ") & (!" + bqmask + ")"
                QM_atoms = len(list(self.traj[mask].top.atoms))
                # Count also link atoms NOT in the Bq region
                for cnt, icrd in enumerate(self.linkobj.link_pos):
                    qm_id  = self.linkobj.qmmm_bonds[cnt][0]
                    crd_qm = self.traj.xyz[igeom][qm_id] 
                    # Assess whether the QM partner is in the Bq region
                    if(list(crd_qm) not in bqcrd):
                        QM_atoms += 1
                    else:
                        pass

                # Now print xyz for atoms NOT in the Bq region
                f.write("{:d}\n\n".format(QM_atoms))

                for cnt, icrd in enumerate(self.traj[self.qmmask].xyz[igeom]):
                    if(list(icrd) in bqcrd):
                        # Iterate, but do NOT print
                        iatom = ATOMNUM[next(itatoms).atomic_number].upper()
                    else:
                    # Print atoms not in bqmask only
                        iatom = ATOMNUM[next(itatoms).atomic_number].upper()
                        f.write(fmt.format(iatom, *icrd))
                # Now print eventual Link atoms, accounting for
                # those NOT present in the Bq region
                for cnt, icrd in enumerate(self.linkobj.link_pos):
                    qm_id  = self.linkobj.qmmm_bonds[cnt][0]
                    crd_qm = self.traj.xyz[igeom][qm_id] 
                    # Assess whether the QM partner is in the Bq region
                    if(list(crd_qm) not in bqcrd):
                        link_at = self.linkobj.link_atom
                        f.write(fmt.format(link_at, *icrd))
#                        link_at = self.linkobj.link_atom + "-Bq"
                    else:
                        pass

                f.write("\n")
        # No Bq atoms present
        else:
            outfile = "geom{:d}.xyz".format(igeom)
            with open(outfile, "w") as f:
                qmatoms = self.traj[self.qmmask].xyz[igeom]
                # Count QM atoms including link atoms
                QM_atoms = self.linkobj.N_qm + self.linkobj.N_link

                f.write("{:d}\n\n".format(QM_atoms))
                for cnt, icrd in enumerate(qmatoms):
                    iatom = ATOMNUM[next(itatoms).atomic_number].upper()
                    f.write(fmt.format(iatom, *icrd))
                # Now print eventual Link atoms
                for cnt, icrd in enumerate(self.linkobj.link_pos):
                    qm_id   = self.linkobj.qmmm_bonds[cnt][0]
                    crd_qm  = self.traj.xyz[igeom][qm_id] 
                    link_at = self.linkobj.link_atom
                    f.write(fmt.format(link_at, *icrd))
                f.write("\n")

    def _write_point_charges(self, f, igeom, prename = None, null_charges = None):
        """Write external point charges
        print_chg = whether to generate a point charges file
        null_charges     = if True, set to zero the point charges"""
        
#        try:
#            print_chg = (self.namelist["externchg"][0] == "true")
#        except:
#            print("No point charges evidenced. Do not generate a point charge file")
#
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
                fmt = "{:<15.6f}{:<15.6f}{:<15.6f}{:<15.6f}{:<10.1f}{:<10.1f}{:<10.1f}\n"
                if(prename != None):
                    chgfile = "charges_{:s}_geom{:d}.xyz".format(prename,igeom)
                else:
                    chgfile = "charges_geom{:d}.xyz".format(igeom)

                # Write charges file
                with open(chgfile, "w") as g:
                    # Count MM atoms accounting for eventual link atoms
                    MM_atoms = self.linkobj.N_mm - self.linkobj.N_link
                    g.write("{:d}\n".format(MM_atoms))
#                    for cnt, icrd in enumerate(self.traj[self.chgmask].xyz[igeom]):
##                        ichg = self.top[self.chgmask].charge[cnt]
#                        ichg = charges[cnt]
#                        g.write(fmt.format(*icrd, ichg, 0.0, 0.0, 0.0))
                    for nat, iat in enumerate(self.traj.top.atoms):
                        if(self.linkobj.labels[nat] == 1):
                            # Consider only MM atoms
                            if(nat not in self.linkobj.all_link):
                                # MM atoms not involved in QM/MM boundary
                                crd = self.traj.xyz[igeom][nat]
                                chg = charges[nat]
                                g.write(fmt.format(*crd, chg, 0.0, 0.0, 0.0))
                            elif(nat in self.linkobj.qmmm_bonds[:,1]):
                                # Print NNs for each QM/MM atom with smeared charges
                                chg_parz = charges[nat]/\
                                           len(self.linkobj.nn_atoms[nat])
                                
                                for inn in self.linkobj.nn_atoms[nat]:
                                    # Ith nearest neighbor
                                    crd = self.traj.xyz[igeom][inn]
                                    chg = charges[inn] + chg_parz
                                    g.write(fmt.format(*crd, chg, 0.0, 0.0, 0.0))
                            else:
                                pass
        else:
            pass

