import numpy as np
import pytraj as pt
from copy import deepcopy
from collections import OrderedDict
from argparse import ArgumentParser
import os

class QMtemplate(object):

    def __init__(self):
        """Parser of the QM template file. General parser"""
        self.namelist = { "header": None,
                          "chgspin": None,
                          "fragments": None,
                          "basis": None,
                          "bsse": None,
                          "externchg": None
                         }

    def parse(self, infile):
        with open(infile, "r") as f:
            lines  = f.read().splitlines()
        nspace = [lines[cnt].replace("&", "").strip().lower()\
            for cnt, line in enumerate(lines) \
            if ("&" in line) and ("end" not in line)]
        
        # args are the general arguments for each template section
        # Parse arguments depending on whether we are defining either 
        # a parent or a child class
        
        local_args = list(self.namelist.keys())
        # 17-01-2022: Set as None attributes all local args
        for iloc in local_args:
            self[iloc] = self.namelist[iloc]
        args       = np.intersect1d(local_args, nspace)

        for name in args:
            self.namelist[name] = self.parse_section(name, lines)
            self[name] = self.parse_section(name, lines)

    def parse_section(self, name, lines):
        """Parse a specific section of the template input.
           Return its contents as a list of strings        
        """

        uname = "&" + name # e.g: &dft
        section = []
#        print("Parsing section: " + uname)
        for line in lines[lines.index(uname) + 1:]:
            row = line.split()
            if(("&" in line) or (len(row) == 0)):
                break;
            else:
                section.append(line)
        return(section)

    def get_fragments(self):
        """Get the fragment mask for each monomer"""
        has_fragments = self.fragments != None
        has_bsse      = self.bsse != None
        has_solvent   = self.solvent != None
        has_closest   = self.closest != None
        if(has_fragments):
            # Get each fragment as a numbered attribute called ._monN
            # where N = 1,2,...
            monomers = []
            for cnt, ifr in enumerate(self.fragments):
                monomers.append(ifr[ifr.index("=")+1:].replace(" ", ""))
#                self["_mon" + cnt] = ifr[ifr.index("=")+1:].replace(" ", "")
        elif(has_bsse):

            for cnt, ifr in enumerate(self.bsse):
                monomers.append(ifr[ifr.index("=")+1:].replace(" ", ""))
#                self["_mon" + cnt] = ifr[ifr.index("=")+1:].replace(" ", "")
        self["monomers"] = monomers

    def __setattr__(self, name, value):
        super(QMtemplate, self).__setattr__(name, value)

    def __setitem__(self, name, value):
        self.__setattr__(name, value)

    def __getattribute__(self, name):
        return(object.__getattribute__(self, name))

    def __getitem__(self, name):
        return(self.__getattribute__(name))

# Decorator for the write_input method(s)
def write_input_bsse(funct):
    def wrapper(*args, **kwargs):
        # args[0] corresponds to self
#        print("args[0] = ", args[0])
        bsse    = args[0].bsse
        chgspin = args[0].chgspin

        if(bsse != None):
            mons  = ["mon1", "mon2"]
            bqidx = [1, 0] # if mon1 -> bqmask = mon2 and vice versa
            null_charges = [True, False]
            # Call write_input method for monomers
            for cnt, imon in enumerate(mons):
                idc = cnt + 1
                charge = chgspin[0].split(",")[2*idc]
                spin   = chgspin[0].split(",")[2*idc + 1]
                bqmask = bsse[bqidx[cnt]][bsse[bqidx[cnt]].index("=")+1:].replace(" ", "")
                print("bqmask ", imon, bqmask)
                funct(*args, **kwargs, 
                      charge  = charge, 
                      spin    = spin,
                      null_charges = null_charges[cnt],
                      prename = imon,
                      bqmask = bqmask) 
            # Call write_input method for complex
            charge = chgspin[0].split(",")[0]
            spin   = chgspin[0].split(",")[1]
            funct(*args, **kwargs,
                  charge  = charge, 
                  spin    = spin,
                  prename = "complex")
        else:
            charge = chgspin[0].split(",")[0]
            spin   = chgspin[0].split(",")[1]
            funct(*args, **kwargs,
                  charge  = charge, 
                  spin    = spin)
    return(wrapper)

