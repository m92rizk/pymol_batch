#!/usr/bin/python

import os, sys, re, shutil, pathlib
from pathlib import Path
import copy
from optparse import OptionParser

dir = os.getcwd()

parser = OptionParser()

parser.add_option("-l", "--ligand", dest="lig")
parser.add_option("-c", "--chain", dest="ch", default="A")
parser.add_option("-p", "--model", dest="pdb")
parser.add_option("", "--apo", action="store_true", dest="apo", default=False)
parser.add_option("", "--compound", action="store_true", dest="compound", default=False)



(options, args) = parser.parse_args()

if options.lig is None:
    lig = raw_input ("Enter ligand (or amino acid) name: ")
else:
    lig = options.lig

LIG = lig.upper()
apo = options.apo
compound = options.compound

if options.pdb is None:
    pdb = raw_input ("Enter pdb name: ")
else:
    pdb = options.pdb
pdbname=Path(pdb).stem

if os.path.isfile('./apo_'+str(pdbname)+'.pdb') and apo is True :
   os.remove('./apo_'+str(pdbname)+'.pdb')

if os.path.isfile(str(LIG)+'.pdb') and compound is True :
   os.remove(str(LIG)+'.pdb')  
    
file = open(pdb, "r")
length = 0
for line in file:
    if (re.search(str(LIG), line) and re.search("HETATM", line)) or (re.search(str(LIG), line) and re.search("ATOM", line)):
        val = str(line[21:22]) # val is the chain
        if str(val) != str(options.ch) and apo is True:
           with open('apo_'+str(pdbname)+'.pdb', "a") as outfile:
               outfile.write(str(line))
        elif str(val) == str(options.ch) and compound is True:
            with open(str(LIG)+'.pdb', "a") as outfile:
               outfile.write(str(line))
    else:
        if apo is True:
            with open('apo_'+str(pdbname)+'.pdb', "a") as outfile:
                outfile.write(str(line))
    if not (re.search("HETATM", line) or re.search("ATOM", line)):
        if compound is True:
            with open(str(LIG)+'.pdb', "a") as outfile:
               outfile.write(str(line))
