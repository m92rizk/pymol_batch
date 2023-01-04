#!/usr/bin/python
# rizk@esrf.fr

import os, sys, re, shutil
import numpy as np
import copy
from optparse import OptionParser


parser = OptionParser()
parser.add_option("-d", "--directory", dest="directory", default=".", help="by default it takes current current directory")
parser.add_option("-f", "--filename", dest="hklfile", default="merged.hkl", help="looking for merged.hkl file by default")
parser.add_option("-u", "--unit_cell", dest="unit", default="", help="unit cell parameters, by default it will be guessed by the program")
parser.add_option("-s", "--space_group", dest="sg", default="", help="spacegroup, by default it will be guessed by the program")
(options, args) = parser.parse_args()


dir = os.getcwd()

if options.unit == "" :
        unit=""
else:
        unit="UNIT_CELL_CONSTANTS=    "+str(options.unit)+"\n"

if options.sg == "" :
        sg=""
else:
        sg="SPACE_GROUP_NUMBER=    "+str(options.sg)+"\n"


if options.directory == "." :
        if os.path.isfile(str(options.hklfile)):
                hklfile = str(options.hklfile)
                print("Found "+str(hklfile))
        else:
                print("Couldn't find any hkl file")

        with open('XDSCONV.INP', 'w') as outxc:
                outxc.write("INPUT_FILE="+str(hklfile)+"\n"+str(unit)+str(sg)+"INCLUDE_RESOLUTION_RANGE=50 1\nOUTPUT_FILE=temp.hkl  CCP4\nFRIEDEL'S_LAW=FALSE\nGENERATE_FRACTION_OF_TEST_REFLECTIONS=0.05")
        os.system('xdsconv >/dev/null')
        os.system('f2mtz HKLOUT temp.mtz<F2MTZ.INP >/dev/null')
        os.system('cad HKLIN1 temp.mtz HKLOUT junk_xdsconv.mtz<<EOF >/dev/null\nLABIN FILE 1 E1=FP E2=SIGFP E3=DANO E4=SIGDANO E5=ISYM E6=FreeRflag\nLABOUT FILE 1 E1=FP E2=SIGFP E3=DANO_sulf E4=SIGDANO_sulf E5=ISYM_sulf E6=FreeRflag\nEND\nEOF ')
        print("done converting hkl to mtz here")
else:
        for folder in os.listdir(dir):
                if folder.startswith(str(options.directory)):
                        os.chdir(str(dir)+'/'+str(folder))
                        if os.path.isfile(str(options.hklfile)):
                                hklfile = str(options.hklfile)
                                print("Found "+str(hklfile))
                        else:
                                print("Couldn't find any hkl file")

                        with open('XDSCONV.INP', 'w') as outxc:
                                outxc.write("INPUT_FILE="+str(hklfile)+"\n"+str(unit)+str(sg)+"INCLUDE_RESOLUTION_RANGE=50 1\nOUTPUT_FILE=temp.hkl  CCP4\nFRIEDEL'S_LAW=FALSE\nGENERATE_FRACTION_OF_TEST_REFLECTIONS=0.05")
                        os.system('xdsconv >/dev/null')
                        os.system('f2mtz HKLOUT temp.mtz<F2MTZ.INP >/dev/null')
                        os.system('cad HKLIN1 temp.mtz HKLOUT junk_xdsconv.mtz<<EOF >/dev/null\nLABIN FILE 1 E1=FP E2=SIGFP E3=DANO E4=SIGDANO E5=ISYM E6=FreeRflag\nLABOUT FILE 1 E1=FP E2=SIGFP E3=DANO_sulf E4=SIGDANO_sulf E5=ISYM_sulf E6=FreeRflag\nEND\nEOF ')
                        print("done converting hkl to mtz in "+str(folder))


print("you may now run get_phases_omit-full.py")
