#!/usr/bin/python3.8

import os, sys, re, shutil
from pathlib import Path
import numpy as np
import time
import copy
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-l", "--ligand", dest="lig", help="ligand name, 3 letters")
parser.add_option("-p", "--model", dest="pdb" , help="pdb file")
parser.add_option("-d", "--directory", dest="directory", default=".", help="leaving this empty, will take current directory")
parser.add_option("-f", "--filename", dest="hklfile", default="merged.hkl", help="looking for merged.hkl file by default")
parser.add_option("", "--phenix", action="store_true", dest="phenix", default=False)
parser.add_option("", "--refmac", action="store_true", dest="refmac", default=False)
parser.add_option("", "--full_map", action="store_true", dest="fullmap", default=False)


(options, args) = parser.parse_args()

if options.lig is None:
    lig = raw_input ("Enter ligand (or amino acid) name: ")
else:
    lig = options.lig

LIG = lig.upper()


if options.pdb is None:
    pdb = raw_input ("Enter pdb name: ")
else:
    pdb = options.pdb

dir = os.getcwd()

if options.fullmap is True:
    omit = ""
    mapname = "full_"
else:
    omit = "apo_"
    mapname = "apo_"

pdbwithpath = pdb
pdb = os.path.basename(pdb)
pdbname = Path(pdb).stem
    
if options.directory == "." :
    if os.path.isfile(str(pdbwithpath)):
        print(str(pdbwithpath)+' loaded ..')
        if not os.path.isfile('../'+str(pdb)):
            shutil.copy(str(pdbwithpath), '../'+str(pdb))
    else:
        print(str(pdbwithpath)+" does not exist\n\nDOWNLOADING !!! "+str(pdb)+"\n\n")
        time.sleep(8)
        os.system('wget http://www.rcsb.org/pdb/files/'+str(pdb))
        os.system('grep -v ANISOU '+str(pdb)+' > renamethis.pdb')
        shutil.move('renamethis.pdb', str(pdb))
        shutil.move(str(pdb), '../'+str(pdb))
    if os.path.isfile('../apo_'+str(pdb)):
        print('apo_'+str(pdb)+' already exists')
    else:
        os.system('./make_apo_lig.py -l '+str(LIG)+' -p ../'+str(pdb)+' --apo')
        shutil.move('apo_'+str(pdb), '../apo_'+str(pdb))
    if os.path.isfile('junk_xdsconv.mtz'):
        print("Found junk_xdsconv.mtz")
    else:
        os.system("./convert_mtz_hkl.py -f "+str(options.hklfile))
    if options.refmac is True:
        os.system('refmac5 hklin junk_xdsconv.mtz xyzin ../'+str(omit)+str(pdb)+' xyzout '+str(omit)+'refmacout_'+str(pdb)+' hklout '+str(mapname)+'map.mtz << eor \nncyc5\nmode rigid\nRIGIdbody NCYCle 5\neor')
    if options.phenix is True:
        os.system('phenix.refine junk_xdsconv.mtz ../'+str(omit)+str(pdb)+'  --overwrite')
        os.system('mv apo_'+str(pdbname)+'_refine_001.mtz '+str(mapname)+'map.mtz')

   
else:
    if os.path.isfile(str(pdbwithpath)):
        print(str(pdbwithpath)+' loaded ..')
        if not os.path.isfile(str(pdb)):
            shutil.copy(str(pdbwithpath), str(pdb))
    else:
        print(str(pdbwithpath)+" does not exist\nDOWNLOADING !!! "+str(pdb)+"\n\n")
        time.sleep(8)
        os.system('wget http://www.rcsb.org/pdb/files/'+str(pdb))
        os.system('grep -v ANISOU '+str(pdb)+' > renamethis.pdb')
        shutil.move('renamethis.pdb', str(pdb))
    if os.path.isfile('apo_'+str(pdb)):
        print('apo_'+str(pdb)+' already exists')
    else:
        os.system('./make_apo_lig.py -l '+str(LIG)+' -p '+str(pdb)+' --apo')
    for folder in os.listdir(dir):
        if folder.startswith(str(options.directory)):
            os.chdir(str(dir)+'/'+str(folder))
            if os.path.isfile('junk_xdsconv.mtz'):
                print("Found junk_xdsconv.mtz")
            else:
                os.system("./convert_mtz_hkl.py -f "+str(options.hklfile))
            if options.refmac is True:
                os.system('refmac5 hklin junk_xdsconv.mtz xyzin ../'+str(omit)+str(pdb)+' xyzout  '+str(omit)+'refmacout_'+str(pdb)+' hklout '+str(mapname)+'map.mtz << eor \nncyc5\nmode rigid\nRIGIdbody NCYCle 5\neor')
            if options.phenix is True:
                os.system('phenix.refine junk_xdsconv.mtz ../'+str(omit)+str(pdb)+'  --overwrite')
                os.system('mv apo_'+str(pdb)+'_refine_001.mtz '+str(mapname)+'map.mtz')
